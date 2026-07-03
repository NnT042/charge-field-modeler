extends Node3D
#
# Atom mode: sim stepping, 3D rendering, scenario spawning.
# UI lives in scenes/ui/atom_hud.tscn (atom_hud.gd) and calls the public
# methods below. All physics is in Rust (AtomSim / atom_core.rs).
#

@export var camera_rig_path: NodePath

var atom_sim: Node  # AtomSim (Rust GDExtension)
var pole_lines: MeshInstance3D
var camera_rig: Node3D
var tuning_panel: CanvasLayer

var profile_ids := {}  # "proton" -> int
var type_renderers := {}  # "proton" -> MultiMeshInstance3D
var substeps_per_frame := 20
var paused := false
var show_force_profile := true
var current_scenario := "protons"

# Force profile visualization
var _fp_mesh: ImmediateMesh
var _fp_mi: MeshInstance3D

# VFX sprinkler
var vfx_mmi: MultiMeshInstance3D
const VFX_EMIT_PER_PARTICLE := 8
const VFX_SPEED := 4.0
const VFX_LIFETIME := 0.4

# Wall-riding orbit geometry (matches Rust atom_scenarios.rs: the M1 orbit
# settles at r=1.3, θ≈11°, v_tan = corotation speed 0.25).
const RIDE_AXIAL := 1.276
const RIDE_LATERAL := 0.248
const COROT_V := 0.25

# ── Lifecycle ──────────────────────────────────────────────────────────

func _ready():
	atom_sim = $AtomSim
	pole_lines = $PoleIndicators
	if camera_rig_path:
		camera_rig = get_node(camera_rig_path)
	_create_tuning_panel()
	_load_profiles()
	_create_envelope_renderers()
	_create_vfx_renderer()
	_create_force_profile_vis()
	spawn_scenario("protons")
	atom_sim.set_running(true)

func _create_tuning_panel() -> void:
	# DEBUG-ONLY: hidden by default (constants are locked and derived — see
	# Couplings::default() in rust/src/atom_core.rs). The HUD's Debug toggle
	# shows it for physics experiments.
	tuning_panel = preload("res://scripts/atom_tuning_panel.gd").new()
	add_child(tuning_panel)
	tuning_panel.setup(atom_sim)
	tuning_panel.visible = false

func set_tuning_visible(on: bool) -> void:
	tuning_panel.visible = on

func _load_profiles():
	var proton = _read_csv("res://config/histogram_proton.csv")
	profile_ids["proton"] = atom_sim.register_profile("proton", 1.0, 1.0, proton)

	var neutron = _read_csv("res://config/histogram_neutron.csv")
	profile_ids["neutron"] = atom_sim.register_profile("neutron", 1.0, 1.0, neutron)

	var electron = _read_csv("res://config/histogram_electron.csv")
	profile_ids["electron"] = atom_sim.register_profile("electron", 1.0 / 1836.0, 0.3, electron)

func _create_envelope_renderers():
	var envelope_mat := StandardMaterial3D.new()
	envelope_mat.vertex_color_use_as_albedo = true
	envelope_mat.transparency = BaseMaterial3D.TRANSPARENCY_ALPHA
	envelope_mat.cull_mode = BaseMaterial3D.CULL_DISABLED
	envelope_mat.shading_mode = BaseMaterial3D.SHADING_MODE_UNSHADED

	for type_name in profile_ids:
		var pid: int = profile_ids[type_name]
		var mesh_buf: PackedFloat32Array = atom_sim.build_profile_mesh(pid, 32, 24)
		var mesh := _build_array_mesh(mesh_buf, envelope_mat)

		var mm := MultiMesh.new()
		mm.transform_format = MultiMesh.TRANSFORM_3D
		mm.use_colors = true
		mm.mesh = mesh

		var mmi := MultiMeshInstance3D.new()
		mmi.multimesh = mm
		add_child(mmi)
		type_renderers[type_name] = mmi

func _build_array_mesh(buf: PackedFloat32Array, mat: Material) -> ArrayMesh:
	if buf.size() < 2:
		return ArrayMesh.new()
	var vert_count := int(buf[0])
	var idx_count := int(buf[1])

	var vertices := PackedVector3Array()
	var normals := PackedVector3Array()
	var colors := PackedColorArray()
	vertices.resize(vert_count)
	normals.resize(vert_count)
	colors.resize(vert_count)

	var offset := 2
	for i in range(vert_count):
		var b := offset + i * 10
		vertices[i] = Vector3(buf[b], buf[b + 1], buf[b + 2])
		normals[i] = Vector3(buf[b + 3], buf[b + 4], buf[b + 5])
		colors[i] = Color(buf[b + 6], buf[b + 7], buf[b + 8], buf[b + 9])

	offset += vert_count * 10
	var indices := PackedInt32Array()
	indices.resize(idx_count)
	for i in range(idx_count):
		indices[i] = int(buf[offset + i])

	var arrays := []
	arrays.resize(Mesh.ARRAY_MAX)
	arrays[Mesh.ARRAY_VERTEX] = vertices
	arrays[Mesh.ARRAY_NORMAL] = normals
	arrays[Mesh.ARRAY_COLOR] = colors
	arrays[Mesh.ARRAY_INDEX] = indices

	var mesh := ArrayMesh.new()
	mesh.add_surface_from_arrays(Mesh.PRIMITIVE_TRIANGLES, arrays)
	mesh.surface_set_material(0, mat)
	return mesh

# ── Scenarios (public API for the HUD) ─────────────────────────────────

func toggle_pause() -> void:
	paused = not paused

func reset_scenario() -> void:
	spawn_scenario(current_scenario)

func spawn_scenario(scenario: String) -> void:
	current_scenario = scenario
	atom_sim.clear_particles()
	match scenario:
		"protons":
			# Two free protons, poles parallel: the stream cushion stands
			# them off — no fusion at ambient pressure.
			atom_sim.spawn_particle(profile_ids["proton"], Vector3(-3, 0, 0), Vector3.ZERO, Vector3.UP)
			atom_sim.spawn_particle(profile_ids["proton"], Vector3(3, 0, 0), Vector3.ZERO, Vector3.UP)
		"hydrogen":
			# Free electron released over the pole: capture, then the
			# wall-riding "circling the drain" orbit.
			atom_sim.spawn_particle(profile_ids["proton"], Vector3.ZERO, Vector3.ZERO, Vector3.UP)
			atom_sim.spawn_particle(profile_ids["electron"],
				Vector3(0.5, 2.0, 0), Vector3(0, 0, 0.5), Vector3.UP)
		"h2_bond":
			# Two H stacked on the pole axis, electrons on the OUTER poles:
			# bare facing poles → stream-cushion bond with vibration.
			_spawn_formed_h(Vector3(0, -3, 0), -1.0)
			_spawn_formed_h(Vector3(0, 3, 0), 1.0)
		"h2_repel":
			# Electrons BETWEEN the protons: stoppered channels, no bond —
			# the pair is driven out past molecular range.
			_spawn_formed_h(Vector3(0, -3, 0), 1.0)
			_spawn_formed_h(Vector3(0, 3, 0), -1.0)
		"alpha", "carbon", "nitrogen", "oxygen":
			atom_sim.spawn_preset(scenario, Vector3.ZERO, Vector3.ZERO, Vector3.UP)
		"helium":
			atom_sim.spawn_preset("alpha", Vector3.ZERO, Vector3.ZERO, Vector3.UP)
			# Electrons riding the two outer proton poles (protons at y=±0.9).
			atom_sim.spawn_particle(profile_ids["electron"],
				Vector3(RIDE_LATERAL, 0.9 + RIDE_AXIAL, 0),
				Vector3(0, 0, -COROT_V), Vector3.UP)
			atom_sim.spawn_particle(profile_ids["electron"],
				Vector3(RIDE_LATERAL, -0.9 - RIDE_AXIAL, 0),
				Vector3(0, 0, -COROT_V), Vector3.DOWN)
		_:
			push_warning("atom_mode: unknown scenario " + scenario)

## One hydrogen atom with the electron already wall-riding.
## `electron_side` +1 = electron on the +Y pole, −1 = on the −Y pole.
func _spawn_formed_h(p_pos: Vector3, electron_side: float) -> void:
	atom_sim.spawn_particle(profile_ids["proton"], p_pos, Vector3.ZERO, Vector3.UP)
	var e_pos := p_pos + Vector3(RIDE_LATERAL, electron_side * RIDE_AXIAL, 0)
	var e_pole := Vector3(0, electron_side, 0)
	atom_sim.spawn_particle(profile_ids["electron"], e_pos, Vector3(0, 0, -COROT_V), e_pole)

func spawn_free(type_name: String) -> void:
	var pid = profile_ids.get(type_name, -1)
	if pid < 0:
		return
	var pos = Vector3(randf_range(-2, 2), randf_range(-2, 2), randf_range(-2, 2))
	atom_sim.spawn_particle(pid, pos, Vector3.ZERO, Vector3.UP)

func _read_csv(path: String) -> PackedFloat32Array:
	var file = FileAccess.open(path, FileAccess.READ)
	if not file:
		push_error("atom_mode: failed to open " + path)
		return PackedFloat32Array()
	var data := PackedFloat32Array()
	var is_header := true
	while not file.eof_reached():
		var line := file.get_line().strip_edges()
		if line.is_empty():
			continue
		if is_header:
			is_header = false
			continue
		var parts := line.split(",")
		if parts.size() >= 2:
			data.append(float(parts[1]))
	return data

func _create_vfx_renderer() -> void:
	var mat := StandardMaterial3D.new()
	mat.shading_mode = BaseMaterial3D.SHADING_MODE_UNSHADED
	mat.vertex_color_use_as_albedo = true
	mat.transparency = BaseMaterial3D.TRANSPARENCY_ALPHA
	mat.billboard_mode = BaseMaterial3D.BILLBOARD_ENABLED
	mat.no_depth_test = true

	var mesh := QuadMesh.new()
	mesh.size = Vector2(0.06, 0.06)
	mesh.material = mat

	var mm := MultiMesh.new()
	mm.transform_format = MultiMesh.TRANSFORM_3D
	mm.use_colors = true
	mm.mesh = mesh

	vfx_mmi = MultiMeshInstance3D.new()
	vfx_mmi.multimesh = mm
	add_child(vfx_mmi)

func _create_force_profile_vis() -> void:
	var mat := StandardMaterial3D.new()
	mat.shading_mode = BaseMaterial3D.SHADING_MODE_UNSHADED
	mat.vertex_color_use_as_albedo = true
	mat.transparency = BaseMaterial3D.TRANSPARENCY_ALPHA
	mat.no_depth_test = true
	_fp_mesh = ImmediateMesh.new()
	_fp_mi = MeshInstance3D.new()
	_fp_mi.mesh = _fp_mesh
	_fp_mi.material_override = mat
	add_child(_fp_mi)

func _draw_force_profiles() -> void:
	_fp_mesh.clear_surfaces()
	if not show_force_profile:
		return
	var n: int = atom_sim.get_particle_count()
	var base_r: float = 0.2
	var vis_base: float = 2.0
	var pts: int = 72
	var pids: PackedInt32Array = atom_sim.get_profile_ids()
	for i in range(n):
		var mass: float = atom_sim.get_profile_mass(pids[i])
		var mass_scale: float = maxf(mass ** 0.25, 0.08)
		var sc: float = vis_base * mass_scale
		for plane_idx in range(2):
			# Emission profile — red/orange (linear, power=1)
			var em: PackedVector3Array = atom_sim.get_profile_ring(
				i, false, base_r, sc, pts, plane_idx, 1.0)
			if em.size() > 1:
				_fp_mesh.surface_begin(Mesh.PRIMITIVE_LINE_STRIP)
				_fp_mesh.surface_set_color(Color(1.0, 0.35, 0.15, 0.55))
				for v in em:
					_fp_mesh.surface_add_vertex(v)
				_fp_mesh.surface_end()
			# Intake profile — blue/cyan (squared, power=2 to match force)
			var ab: PackedVector3Array = atom_sim.get_profile_ring(
				i, true, base_r, sc, pts, plane_idx, 2.0)
			if ab.size() > 1:
				_fp_mesh.surface_begin(Mesh.PRIMITIVE_LINE_STRIP)
				_fp_mesh.surface_set_color(Color(0.2, 0.55, 1.0, 0.55))
				for v in ab:
					_fp_mesh.surface_add_vertex(v)
				_fp_mesh.surface_end()

# ── Frame loop ─────────────────────────────────────────────────────────

func _process(_delta):
	if not paused and atom_sim.is_running():
		atom_sim.step_n(substeps_per_frame)
	_update_rendering()
	_update_vfx(_delta)
	_draw_force_profiles()

func _update_vfx(delta: float) -> void:
	if paused or not atom_sim.is_running() or not atom_sim.is_vfx_enabled():
		vfx_mmi.multimesh.instance_count = 0
		return
	var buf: PackedFloat32Array = atom_sim.advance_vfx(
		delta, VFX_EMIT_PER_PARTICLE, VFX_SPEED, VFX_LIFETIME
	)
	var count: int = buf.size() / 16
	vfx_mmi.multimesh.instance_count = count
	if count > 0:
		vfx_mmi.multimesh.set_buffer(buf)

# ── Input ──────────────────────────────────────────────────────────────

func _input(event: InputEvent) -> void:
	if not (event is InputEventKey) or not event.pressed or event.echo:
		return
	if (event as InputEventKey).keycode == KEY_TAB:
		get_viewport().set_input_as_handled()
		get_tree().change_scene_to_file("res://scenes/main.tscn")

func _unhandled_key_input(event: InputEvent):
	if not (event is InputEventKey):
		return
	var k: InputEventKey = event
	if not k.pressed or k.echo:
		return
	var handled := true
	match k.keycode:
		KEY_SPACE:
			toggle_pause()
		KEY_BACKSPACE:
			reset_scenario()
		KEY_P:
			spawn_free("proton")
		KEY_N:
			spawn_free("neutron")
		KEY_E:
			spawn_free("electron")
		KEY_H:
			spawn_scenario("hydrogen")
		KEY_D:
			spawn_scenario("h2_bond")
		KEY_F:
			show_force_profile = not show_force_profile
		KEY_V:
			atom_sim.set_vfx_enabled(not atom_sim.is_vfx_enabled())
		KEY_BRACKETRIGHT:
			substeps_per_frame = mini(substeps_per_frame + 5, 200)
		KEY_BRACKETLEFT:
			substeps_per_frame = maxi(substeps_per_frame - 5, 1)
		_:
			handled = false
	if handled:
		get_viewport().set_input_as_handled()

# ── Rendering ──────────────────────────────────────────────────────────

func _update_rendering():
	# Per-type envelope mesh rendering
	for type_name in type_renderers:
		var pid: int = profile_ids[type_name]
		var count: int = atom_sim.count_particles_with_profile(pid)
		var mmi: MultiMeshInstance3D = type_renderers[type_name]
		if count == 0:
			mmi.multimesh.instance_count = 0
			continue
		if mmi.multimesh.instance_count != count:
			mmi.multimesh.instance_count = count
		var buf: PackedFloat32Array = atom_sim.build_multimesh_buffer_for_profile(pid)
		mmi.multimesh.set_buffer(buf)

	# Pole indicator lines
	var total: int = atom_sim.get_particle_count()
	if total == 0:
		pole_lines.mesh = null
		return
	var pole_buf: PackedFloat32Array = atom_sim.build_pole_indicator_buffer()
	_update_pole_mesh(pole_buf, total)

func _update_pole_mesh(buf: PackedFloat32Array, count: int):
	if count == 0:
		pole_lines.mesh = null
		return
	var arrays := []
	arrays.resize(Mesh.ARRAY_MAX)

	var verts := PackedVector3Array()
	var colors := PackedColorArray()
	verts.resize(count * 2)
	colors.resize(count * 2)

	for i in range(count * 2):
		var base := i * 6
		verts[i] = Vector3(buf[base], buf[base + 1], buf[base + 2])
		colors[i] = Color(buf[base + 3], buf[base + 4], buf[base + 5])

	arrays[Mesh.ARRAY_VERTEX] = verts
	arrays[Mesh.ARRAY_COLOR] = colors

	var mesh := ArrayMesh.new()
	mesh.add_surface_from_arrays(Mesh.PRIMITIVE_LINES, arrays)
	pole_lines.mesh = mesh
