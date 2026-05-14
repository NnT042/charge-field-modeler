extends Node3D
## GPU compute orchestrator for the charge field simulation.
## Three-pass pipeline: field_update → collision_detect → impulse_sum.

@export var focus_particle_path: NodePath
@export var field_renderer_path: NodePath
@export var default_photon_count: int = 50000

var _focus: Node3D
var _renderer: Node3D
var _rust_sim: Node

var _rd: RenderingDevice

var _update_shader: RID
var _collide_shader: RID
var _reduce_shader: RID
var _trace_shader: RID
var _update_pipeline: RID
var _collide_pipeline: RID
var _reduce_pipeline: RID
var _trace_pipeline: RID

var _photon_buffer: RID
var _params_buffer: RID
var _focus_state_buffer: RID
var _impulse_buffer: RID
var _reduction_buffer: RID
var _reduce_params_buffer: RID

var _trace_segments_buffer: RID

var _update_set: RID
var _collide_set: RID
var _trace_set: RID
var _reduce_set: RID
var _buffers_valid: bool = false

var _photon_count: int = 0
var _workgroups: int = 0
var _frame_number: int = 0
var _running: bool = false
var _net_impulse: Vector3 = Vector3.ZERO
var _collision_count: int = 0
var _show_exit_holes: bool = false
var _show_collision_pointers: bool = false
var _show_heatmap: bool = false
var _cloud_blanked: bool = false
var _stream_mode: bool = false

var _hole_multimesh: MultiMesh
var _hole_instance: MultiMeshInstance3D

var _pointer_multimesh: MultiMesh
var _pointer_instance: MultiMeshInstance3D

var _heatmap_sphere: MeshInstance3D
var _heatmap_image: Image
var _heatmap_texture: ImageTexture
var _heatmap_material: ShaderMaterial


func _ready() -> void:
	_focus = get_node_or_null(focus_particle_path) as Node3D
	_renderer = get_node_or_null(field_renderer_path) as Node3D

	_rust_sim = FieldSim.new()
	add_child(_rust_sim)

	if _focus:
		var eff_r: float = _focus.call("effective_radius")
		_rust_sim.call("adapt_sim_radius", eff_r)

	_setup_exit_hole_renderer()
	_setup_pointer_renderer()
	_setup_heatmap_sphere()

	_rd = RenderingServer.create_local_rendering_device()
	_update_shader = _load_shader("res://shaders/compute/field_update.glsl")
	_collide_shader = _load_shader("res://shaders/compute/collision_detect.glsl")
	_trace_shader = _load_shader("res://shaders/compute/trace_collide.glsl")
	_reduce_shader = _load_shader("res://shaders/compute/impulse_sum.glsl")
	if not _update_shader.is_valid():
		return
	_update_pipeline = _rd.compute_pipeline_create(_update_shader)
	_collide_pipeline = _rd.compute_pipeline_create(_collide_shader)
	_trace_pipeline = _rd.compute_pipeline_create(_trace_shader)
	_reduce_pipeline = _rd.compute_pipeline_create(_reduce_shader)


func _load_shader(path: String) -> RID:
	var shader_file: RDShaderFile = load(path)
	if shader_file == null:
		push_error("Could not load shader: " + path)
		return RID()
	var spirv: RDShaderSPIRV = shader_file.get_spirv()
	var err: String = spirv.get_stage_compile_error(RenderingDevice.SHADER_STAGE_COMPUTE)
	if err != "":
		push_error(path + " compile error: " + err)
		return RID()
	return _rd.shader_create_from_spirv(spirv)


func _physics_process(delta: float) -> void:
	var time_scale: float = TAU
	if _focus:
		time_scale = _focus.call("get_time_scale")

	# Keep effective_radius synced for disc spawn
	if _focus:
		var eff_r_update: float = _focus.call("effective_radius")
		_rust_sim.call("adapt_sim_radius", eff_r_update)

	# GPU cloud passes (only when cloud is active)
	if _running and _photon_count > 0:
		_frame_number += 1

		var seg_count: int = 0
		if _focus:
			seg_count = int(_focus.call("get_trace_segment_count", 4096))
			if seg_count > 0:
				var trace_buf: PackedFloat32Array = _focus.call("pack_trace_segments", 4096)
				var trace_bytes: PackedByteArray = trace_buf.to_byte_array()
				_rd.buffer_update(_trace_segments_buffer, 0, trace_bytes.size(), trace_bytes)

		var center: Vector3 = Vector3.ZERO
		var params: PackedFloat32Array = _rust_sim.call(
			"pack_sim_params", center, delta, time_scale, _frame_number, seg_count
		)
		_rd.buffer_update(_params_buffer, 0, params.to_byte_array().size(), params.to_byte_array())

		var base_solid: PackedFloat32Array = _focus.call("pack_base_solid") if _focus else PackedFloat32Array([0.0, 0.0, 0.0, 0.0])
		var bs_bytes: PackedByteArray = base_solid.to_byte_array()
		_rd.buffer_update(_focus_state_buffer, 0, bs_bytes.size(), bs_bytes)

		var cl: int = _rd.compute_list_begin()
		_rd.compute_list_bind_compute_pipeline(cl, _update_pipeline)
		_rd.compute_list_bind_uniform_set(cl, _update_set, 0)
		_rd.compute_list_dispatch(cl, _workgroups, 1, 1)
		_rd.compute_list_end()

		cl = _rd.compute_list_begin()
		_rd.compute_list_bind_compute_pipeline(cl, _trace_pipeline)
		_rd.compute_list_bind_uniform_set(cl, _trace_set, 0)
		_rd.compute_list_dispatch(cl, _workgroups, 1, 1)
		_rd.compute_list_end()

		cl = _rd.compute_list_begin()
		_rd.compute_list_bind_compute_pipeline(cl, _reduce_pipeline)
		_rd.compute_list_bind_uniform_set(cl, _reduce_set, 0)
		_rd.compute_list_dispatch(cl, _workgroups, 1, 1)
		_rd.compute_list_end()

		_rd.submit()
		_rd.sync()

		var output: PackedByteArray = _rd.buffer_get_data(_photon_buffer)
		_rust_sim.call("parse_readback", output)

		if _show_exit_holes or _show_collision_pointers or _show_heatmap:
			var impulse_raw: PackedByteArray = _rd.buffer_get_data(_impulse_buffer)
			_rust_sim.call("detect_exit_holes_from_impulses", impulse_raw)

		var reduce_data: PackedByteArray = _rd.buffer_get_data(_reduction_buffer)
		_net_impulse += _rust_sim.call("sum_partial_impulses", reduce_data)
		_collision_count += int(_rust_sim.call("get_last_collision_count"))

	var eff_r: float = 1.0
	if _focus:
		eff_r = float(_focus.call("effective_radius"))

	# Render cloud (only rebuild buffer when sim is running and not blanked)
	if _running and not _cloud_blanked:
		var total: int = int(_rust_sim.call("get_total_display_count"))
		if _renderer and total > 0:
			_renderer.call("set_instance_count", total)
			var mm_buf: PackedFloat32Array = _rust_sim.call("build_multimesh_buffer")
			_renderer.call("update_from_buffer", mm_buf)
		elif _renderer:
			_renderer.call("set_instance_count", 0)

	# Render exit holes on ghost sphere (only rebuild while sim produces new data)
	if _show_exit_holes and _running:
		var hole_count: int = int(_rust_sim.call("get_exit_hole_count"))
		if hole_count > 0:
			_hole_multimesh.instance_count = hole_count
			var hole_buf: PackedFloat32Array = _rust_sim.call("build_exit_holes_buffer", eff_r)
			_hole_multimesh.set_buffer(hole_buf)
		else:
			_hole_multimesh.instance_count = 0
	elif not _show_exit_holes:
		_hole_multimesh.instance_count = 0

	# Render collision pointers (only rebuild while sim produces new data)
	if _show_collision_pointers and _running:
		var ptr_count: int = int(_rust_sim.call("get_collision_pointer_count"))
		if ptr_count > 0:
			_pointer_multimesh.instance_count = ptr_count * 2
			var ptr_buf: PackedFloat32Array = _rust_sim.call("build_collision_pointers_buffer", eff_r)
			_pointer_multimesh.set_buffer(ptr_buf)
		else:
			_pointer_multimesh.instance_count = 0
	elif not _show_collision_pointers:
		_pointer_multimesh.instance_count = 0

	# Update heatmap sphere (every 10 frames to save CPU)
	if _show_heatmap and _running and _frame_number % 10 == 0:
		_update_heatmap_texture(eff_r)


func start_field(count: int = -1) -> void:
	if count > 0:
		_photon_count = count
	elif _photon_count == 0:
		_photon_count = default_photon_count

	if _focus:
		var eff_r: float = _focus.call("effective_radius")
		_rust_sim.call("adapt_sim_radius", eff_r)

	_rust_sim.call("init_field", _photon_count)
	_workgroups = ceili(float(_photon_count) / 256.0)
	_create_buffers()

	if _renderer:
		_renderer.call("set_instance_count", _photon_count)
		var mm_buf: PackedFloat32Array = _rust_sim.call("build_multimesh_buffer")
		_renderer.call("update_from_buffer", mm_buf)

	_running = true
	_frame_number = 0
	if _focus:
		_focus.call("freeze_trace")


func stop_field() -> void:
	_running = false
	if _focus:
		_focus.call("unfreeze_trace")


func reset_field() -> void:
	_running = false
	_photon_count = 0
	_frame_number = 0
	_net_impulse = Vector3.ZERO
	_collision_count = 0
	_cleanup_buffers()
	_rust_sim.call("clear_field")
	if _renderer:
		_renderer.call("set_instance_count", 0)
	if _focus:
		_focus.call("unfreeze_trace")
	_clear_heatmap_visual()


func is_running() -> bool:
	return _running


func get_photon_count() -> int:
	return _photon_count


func get_net_impulse() -> Vector3:
	return _net_impulse


func get_collision_count() -> int:
	return _collision_count


func clear_exit_holes() -> void:
	_rust_sim.call("clear_exit_holes")


func set_photon_ratio(ratio: float) -> void:
	_rust_sim.call("set_photon_ratio", ratio)


func set_field_direction(dir: Vector3) -> void:
	_rust_sim.call("set_field_direction", dir)


func set_direction_strength(strength: float) -> void:
	_rust_sim.call("set_direction_strength", strength)


func set_direction_scatter(scatter: float) -> void:
	_rust_sim.call("set_direction_scatter", scatter)


func set_direction_mask(mask: int) -> void:
	_rust_sim.call("set_direction_mask", mask)


func set_cloud_blanked(blanked: bool) -> void:
	_cloud_blanked = blanked
	if blanked and _renderer:
		_renderer.call("set_visible_field", false)
	elif not blanked and _renderer:
		_renderer.call("set_visible_field", true)


func set_stream_mode(enabled: bool) -> void:
	_stream_mode = enabled
	_rust_sim.call("set_stream_jitter", 1.0 if enabled else 0.0)


func set_show_exit_holes(show: bool) -> void:
	_show_exit_holes = show
	if not show:
		_hole_multimesh.instance_count = 0
		_rust_sim.call("clear_exit_holes")


func set_show_collision_pointers(show: bool) -> void:
	_show_collision_pointers = show
	if not show:
		_pointer_multimesh.instance_count = 0
		_rust_sim.call("clear_collision_pointers")


func set_show_heatmap(show: bool) -> void:
	_show_heatmap = show
	_heatmap_sphere.visible = show
	if not show:
		_rust_sim.call("clear_heatmap")


func get_mean_emission_angle() -> float:
	return float(_rust_sim.call("get_mean_emission_angle"))


func get_heatmap_total_hits() -> int:
	return int(_rust_sim.call("get_heatmap_total_hits"))


func _update_heatmap_texture(eff_r: float) -> void:
	var buf: PackedFloat32Array = _rust_sim.call("get_heatmap_buffer")
	if buf.size() == 0:
		return
	var dims: Vector2 = _rust_sim.call("get_heatmap_dimensions")
	var w: int = int(dims.x)
	var h: int = int(dims.y)

	var byte_buf: PackedByteArray = buf.to_byte_array()
	_heatmap_image = Image.create_from_data(w, h, false, Image.FORMAT_RGBF, byte_buf)
	_heatmap_texture.update(_heatmap_image)

	# Scale sphere to match ghost sphere outer radius
	var scale: float = eff_r * 0.5 + 0.05
	_heatmap_sphere.scale = Vector3(scale, scale, scale)


func _create_buffers() -> void:
	_cleanup_buffers()

	var init_data: PackedByteArray = _rust_sim.call("get_init_buffer")
	_photon_buffer = _rd.storage_buffer_create(init_data.size(), init_data)

	var params_init: PackedByteArray = PackedByteArray()
	params_init.resize(80)
	_params_buffer = _rd.uniform_buffer_create(80, params_init)

	# Focus state: vec4 (16 bytes, but std140 needs 16-byte alignment)
	var focus_init: PackedByteArray = PackedByteArray()
	focus_init.resize(16)
	_focus_state_buffer = _rd.uniform_buffer_create(16, focus_init)

	# Impulse buffer: one vec4 per photon
	var impulse_size: int = _photon_count * 16
	_impulse_buffer = _rd.storage_buffer_create(impulse_size)

	# Reduction buffer: one vec4 per workgroup
	var reduce_size: int = _workgroups * 16
	_reduction_buffer = _rd.storage_buffer_create(reduce_size)

	# Reduction params uniform (element_count + 3 pad)
	var rp: PackedByteArray = PackedByteArray()
	rp.resize(16)
	rp.encode_u32(0, _photon_count)
	_reduce_params_buffer = _rd.uniform_buffer_create(16, rp)

	# Trace segments buffer: 4096 segments * 64 bytes = 256KB
	var trace_size: int = 4096 * 64
	_trace_segments_buffer = _rd.storage_buffer_create(trace_size)

	# --- Uniform sets ---

	# Pass 1: field_update — photons(0) + params(1)
	var u0: RDUniform = RDUniform.new()
	u0.uniform_type = RenderingDevice.UNIFORM_TYPE_STORAGE_BUFFER
	u0.binding = 0
	u0.add_id(_photon_buffer)

	var u1: RDUniform = RDUniform.new()
	u1.uniform_type = RenderingDevice.UNIFORM_TYPE_UNIFORM_BUFFER
	u1.binding = 1
	u1.add_id(_params_buffer)

	_update_set = _rd.uniform_set_create([u0, u1], _update_shader, 0)

	# Pass 2 (legacy): collision_detect — photons(0) + params(1) + focus_state(2) + impulses(3)
	var c0: RDUniform = RDUniform.new()
	c0.uniform_type = RenderingDevice.UNIFORM_TYPE_STORAGE_BUFFER
	c0.binding = 0
	c0.add_id(_photon_buffer)

	var c1: RDUniform = RDUniform.new()
	c1.uniform_type = RenderingDevice.UNIFORM_TYPE_UNIFORM_BUFFER
	c1.binding = 1
	c1.add_id(_params_buffer)

	var c2: RDUniform = RDUniform.new()
	c2.uniform_type = RenderingDevice.UNIFORM_TYPE_UNIFORM_BUFFER
	c2.binding = 2
	c2.add_id(_focus_state_buffer)

	var c3: RDUniform = RDUniform.new()
	c3.uniform_type = RenderingDevice.UNIFORM_TYPE_STORAGE_BUFFER
	c3.binding = 3
	c3.add_id(_impulse_buffer)

	_collide_set = _rd.uniform_set_create([c0, c1, c2, c3], _collide_shader, 0)

	# Pass 2 (trace): trace_collide — photons(0) + params(1) + focus_state(2) + impulses(3) + segments(4)
	var t0: RDUniform = RDUniform.new()
	t0.uniform_type = RenderingDevice.UNIFORM_TYPE_STORAGE_BUFFER
	t0.binding = 0
	t0.add_id(_photon_buffer)

	var t1: RDUniform = RDUniform.new()
	t1.uniform_type = RenderingDevice.UNIFORM_TYPE_UNIFORM_BUFFER
	t1.binding = 1
	t1.add_id(_params_buffer)

	var t2: RDUniform = RDUniform.new()
	t2.uniform_type = RenderingDevice.UNIFORM_TYPE_UNIFORM_BUFFER
	t2.binding = 2
	t2.add_id(_focus_state_buffer)

	var t3: RDUniform = RDUniform.new()
	t3.uniform_type = RenderingDevice.UNIFORM_TYPE_STORAGE_BUFFER
	t3.binding = 3
	t3.add_id(_impulse_buffer)

	var t4: RDUniform = RDUniform.new()
	t4.uniform_type = RenderingDevice.UNIFORM_TYPE_STORAGE_BUFFER
	t4.binding = 4
	t4.add_id(_trace_segments_buffer)

	_trace_set = _rd.uniform_set_create([t0, t1, t2, t3, t4], _trace_shader, 0)

	# Pass 3: impulse_sum — impulses(0) + partial_sums(1) + reduce_params(2)
	var r0: RDUniform = RDUniform.new()
	r0.uniform_type = RenderingDevice.UNIFORM_TYPE_STORAGE_BUFFER
	r0.binding = 0
	r0.add_id(_impulse_buffer)

	var r1: RDUniform = RDUniform.new()
	r1.uniform_type = RenderingDevice.UNIFORM_TYPE_STORAGE_BUFFER
	r1.binding = 1
	r1.add_id(_reduction_buffer)

	var r2: RDUniform = RDUniform.new()
	r2.uniform_type = RenderingDevice.UNIFORM_TYPE_UNIFORM_BUFFER
	r2.binding = 2
	r2.add_id(_reduce_params_buffer)

	_reduce_set = _rd.uniform_set_create([r0, r1, r2], _reduce_shader, 0)
	_buffers_valid = true


func _setup_exit_hole_renderer() -> void:
	var shader: Shader = load("res://shaders/visual/field_dot.gdshader")
	var mat := ShaderMaterial.new()
	mat.shader = shader

	var quad := QuadMesh.new()
	quad.size = Vector2(1.0, 1.0)
	quad.material = mat

	_hole_multimesh = MultiMesh.new()
	_hole_multimesh.transform_format = MultiMesh.TRANSFORM_3D
	_hole_multimesh.use_colors = true
	_hole_multimesh.mesh = quad
	_hole_multimesh.instance_count = 0

	_hole_instance = MultiMeshInstance3D.new()
	_hole_instance.multimesh = _hole_multimesh
	add_child(_hole_instance)


func _setup_pointer_renderer() -> void:
	var shader: Shader = load("res://shaders/visual/field_dot.gdshader")
	var mat := ShaderMaterial.new()
	mat.shader = shader

	var quad := QuadMesh.new()
	quad.size = Vector2(1.0, 1.0)
	quad.material = mat

	_pointer_multimesh = MultiMesh.new()
	_pointer_multimesh.transform_format = MultiMesh.TRANSFORM_3D
	_pointer_multimesh.use_colors = true
	_pointer_multimesh.mesh = quad
	_pointer_multimesh.instance_count = 0

	_pointer_instance = MultiMeshInstance3D.new()
	_pointer_instance.multimesh = _pointer_multimesh
	add_child(_pointer_instance)


func _setup_heatmap_sphere() -> void:
	var shader: Shader = load("res://shaders/visual/heatmap_sphere.gdshader")
	_heatmap_material = ShaderMaterial.new()
	_heatmap_material.shader = shader

	_heatmap_image = Image.create(64, 32, false, Image.FORMAT_RGBF)
	_heatmap_texture = ImageTexture.create_from_image(_heatmap_image)
	_heatmap_material.set_shader_parameter("heatmap_texture", _heatmap_texture)

	var sphere := SphereMesh.new()
	sphere.radius = 1.0
	sphere.height = 2.0
	sphere.radial_segments = 48
	sphere.rings = 24

	_heatmap_sphere = MeshInstance3D.new()
	_heatmap_sphere.mesh = sphere
	_heatmap_sphere.material_override = _heatmap_material
	_heatmap_sphere.cast_shadow = GeometryInstance3D.SHADOW_CASTING_SETTING_OFF
	_heatmap_sphere.visible = false
	add_child(_heatmap_sphere)


func _clear_heatmap_visual() -> void:
	_heatmap_image = Image.create(64, 32, false, Image.FORMAT_RGBF)
	_heatmap_texture.update(_heatmap_image)
	_heatmap_sphere.visible = _show_heatmap


func _cleanup_buffers() -> void:
	if not _buffers_valid:
		return
	_rd.free_rid(_reduce_set)
	_rd.free_rid(_trace_set)
	_rd.free_rid(_collide_set)
	_rd.free_rid(_update_set)
	_rd.free_rid(_reduce_params_buffer)
	_rd.free_rid(_reduction_buffer)
	_rd.free_rid(_impulse_buffer)
	_rd.free_rid(_trace_segments_buffer)
	_rd.free_rid(_focus_state_buffer)
	_rd.free_rid(_params_buffer)
	_rd.free_rid(_photon_buffer)
	_buffers_valid = false


func _exit_tree() -> void:
	_cleanup_buffers()
	if _update_pipeline.is_valid():
		_rd.free_rid(_update_pipeline)
	if _collide_pipeline.is_valid():
		_rd.free_rid(_collide_pipeline)
	if _trace_pipeline.is_valid():
		_rd.free_rid(_trace_pipeline)
	if _reduce_pipeline.is_valid():
		_rd.free_rid(_reduce_pipeline)
	if _update_shader.is_valid():
		_rd.free_rid(_update_shader)
	if _collide_shader.is_valid():
		_rd.free_rid(_collide_shader)
	if _trace_shader.is_valid():
		_rd.free_rid(_trace_shader)
	if _reduce_shader.is_valid():
		_rd.free_rid(_reduce_shader)
	_rd.free()
