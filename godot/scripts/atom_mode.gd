extends Node3D

@export var camera_rig_path: NodePath

var atom_sim: Node  # AtomSim (Rust GDExtension)
var pole_lines: MeshInstance3D
var camera_rig: Node3D
var info_label: Label

var profile_ids := {}  # "proton" -> int
var type_renderers := {}  # "proton" -> MultiMeshInstance3D
var substeps_per_frame := 20
var paused := false

# VFX sprinkler
var vfx_mmi: MultiMeshInstance3D
const VFX_EMIT_PER_PARTICLE := 8
const VFX_SPEED := 4.0
const VFX_LIFETIME := 0.4

# ── Lifecycle ──────────────────────────────────────────────────────────

func _ready():
	atom_sim = $AtomSim
	pole_lines = $PoleIndicators
	if camera_rig_path:
		camera_rig = get_node(camera_rig_path)
	_create_info_label()
	_create_tuning_panel()
	_load_profiles()
	_create_envelope_renderers()
	_create_vfx_renderer()
	_spawn_test_scenario()
	atom_sim.set_running(true)

func _create_info_label():
	var canvas := CanvasLayer.new()
	add_child(canvas)
	info_label = Label.new()
	info_label.position = Vector2(12, 8)
	info_label.add_theme_font_size_override("font_size", 14)
	info_label.add_theme_color_override("font_color", Color(0.8, 0.85, 0.9))
	info_label.add_theme_color_override("font_shadow_color", Color(0, 0, 0, 0.6))
	info_label.add_theme_constant_override("shadow_offset_x", 1)
	info_label.add_theme_constant_override("shadow_offset_y", 1)
	canvas.add_child(info_label)

func _create_tuning_panel() -> void:
	# DISPOSABLE: delete this function + its call in _ready() and the
	# atom_tuning_panel.gd file to remove the slider UI. No physics here.
	var panel = preload("res://scripts/atom_tuning_panel.gd").new()
	add_child(panel)
	panel.setup(atom_sim)

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

func _spawn_test_scenario():
	atom_sim.clear_particles()
	# Two protons, poles up, separated along X
	atom_sim.spawn_particle(profile_ids["proton"], Vector3(-3, 0, 0), Vector3.ZERO, Vector3.UP)
	atom_sim.spawn_particle(profile_ids["proton"], Vector3(3, 0, 0), Vector3.ZERO, Vector3.UP)

func _spawn_hydrogen():
	# Proton at origin, electron near north pole with tangential velocity for orbit
	var p_pos := Vector3(0, 0, 0)
	atom_sim.spawn_particle(profile_ids["proton"], p_pos, Vector3.ZERO, Vector3.UP)
	atom_sim.spawn_particle(profile_ids["electron"],
		p_pos + Vector3(0.5, 2.0, 0),
		Vector3(0, 0, 0.5),
		Vector3.UP)

func _spawn_h2_test():
	# Two hydrogen atoms: protons separated along X, electrons on outer poles
	# Left H: proton at (-4,0,0) pole up, electron at (-4,3,0) above pole
	atom_sim.spawn_particle(profile_ids["proton"], Vector3(-4, 0, 0), Vector3.ZERO, Vector3.UP)
	atom_sim.spawn_particle(profile_ids["electron"], Vector3(-4, 3, 0), Vector3.ZERO, Vector3.UP)
	# Right H: proton at (4,0,0) pole up, electron at (4,3,0) above pole
	atom_sim.spawn_particle(profile_ids["proton"], Vector3(4, 0, 0), Vector3.ZERO, Vector3.UP)
	atom_sim.spawn_particle(profile_ids["electron"], Vector3(4, 3, 0), Vector3.ZERO, Vector3.UP)

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

# ── Frame loop ─────────────────────────────────────────────────────────

func _process(_delta):
	if not paused and atom_sim.is_running():
		atom_sim.step_n(substeps_per_frame)
	_update_rendering()
	_update_vfx(_delta)
	_update_info()

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

func _update_info():
	var n: int = atom_sim.get_particle_count()
	var t: float = atom_sim.get_time()
	var ke: float = atom_sim.get_total_kinetic_energy()
	var status := "PAUSED" if paused else "RUNNING"
	var lines := PackedStringArray()
	lines.append("ATOM MODE [%s]  t=%.4f  KE=%s" % [status, t, String.num_scientific(ke)])
	var dt: float = atom_sim.get_timestep()
	lines.append("Particles: %d   dt=%s   substeps=%d" % [n, String.num_scientific(dt), substeps_per_frame])
	var g_q: float = atom_sim.get_gravity_coupling()
	var c_q: float = atom_sim.get_charge_coupling()
	var p_amb: float = atom_sim.get_ambient_pressure()
	var v_q: float = atom_sim.get_vortex_coupling()
	var d_q: float = atom_sim.get_drag_coupling()
	lines.append("G_q=%.2f  C_q=%.2f  P_amb=%.3f  T_q=%.3f  V_q=%.3f  D_q=%.2f" % [g_q, c_q, p_amb, float(atom_sim.get_torque_coupling()), v_q, d_q])
	var oi: PackedFloat32Array = atom_sim.get_orbit_info()
	if oi.size() >= 5 and oi[0] > 0.0:
		lines.append("Orbit: r=%.2f  v_tan=%.3f  v_rad=%.3f  ->G_q*=%.3f  E_axis=%.4f" % [oi[0], oi[1], oi[2], oi[3], oi[4]])
	for i in range(mini(n, 6)):
		var pos: Vector3 = atom_sim.get_particle_position(i)
		var spd: float = atom_sim.get_particle_speed(i)
		var pname: String = atom_sim.get_profile_name(atom_sim.get_profile_ids()[i])
		lines.append("  #%d %s  pos=(%.2f,%.2f,%.2f)  v=%.3f" % [i, pname, pos.x, pos.y, pos.z, spd])
	if n >= 2:
		lines.append("  d(0,1)=%.4f" % atom_sim.get_pair_distance(0, 1))
	lines.append("")
	var vfx_status := "ON" if atom_sim.is_vfx_enabled() else "OFF"
	lines.append("VFX: %s (%d)  " % [vfx_status, atom_sim.get_vfx_count()])
	lines.append("[Space] pause  [Bksp] reset  [P/N/E] spawn  [H] hydrogen  [D] H2")
	lines.append("[V] VFX  [G/Shift+G] G_q  [C/Shift+C] C_q  [A/Shift+A] P_amb")
	lines.append("[W/Shift+W] V_q  [T/Shift+T] T_q  [R/Shift+R] D_q  [K] calibrate  [Tab] Phase 1")
	info_label.text = "\n".join(lines)

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
	var shift := k.shift_pressed
	match k.keycode:
		KEY_SPACE:
			paused = not paused
		KEY_BACKSPACE:
			atom_sim.clear_particles()
			_spawn_test_scenario()
		KEY_P:
			_spawn_at_cursor("proton")
		KEY_N:
			_spawn_at_cursor("neutron")
		KEY_E:
			_spawn_at_cursor("electron")
		KEY_H:
			atom_sim.clear_particles()
			_spawn_hydrogen()
		KEY_D:
			atom_sim.clear_particles()
			_spawn_h2_test()
		KEY_V:
			atom_sim.set_vfx_enabled(not atom_sim.is_vfx_enabled())
		KEY_G:
			var cur: float = atom_sim.get_gravity_coupling()
			atom_sim.set_gravity_coupling(cur * (0.8 if shift else 1.25))
		KEY_C:
			var cur: float = atom_sim.get_charge_coupling()
			atom_sim.set_charge_coupling(cur * (0.8 if shift else 1.25))
		KEY_A:
			var cur: float = atom_sim.get_ambient_pressure()
			if shift:
				atom_sim.set_ambient_pressure(maxf(cur - 0.01, 0.0))
			else:
				atom_sim.set_ambient_pressure(cur + 0.01)
		KEY_W:
			var cur: float = atom_sim.get_vortex_coupling()
			atom_sim.set_vortex_coupling(cur * (0.8 if shift else 1.25))
		KEY_T:
			var cur: float = atom_sim.get_torque_coupling()
			atom_sim.set_torque_coupling(cur * (0.8 if shift else 1.25))
		KEY_R:
			var cur: float = atom_sim.get_drag_coupling()
			atom_sim.set_drag_coupling(cur * (0.8 if shift else 1.25))
		KEY_K:
			var g: float = atom_sim.auto_calibrate_polar()
			print("[atom] auto_calibrate_polar -> G_q = %.4f" % g)
		KEY_BRACKETRIGHT:
			substeps_per_frame = mini(substeps_per_frame + 5, 200)
		KEY_BRACKETLEFT:
			substeps_per_frame = maxi(substeps_per_frame - 5, 1)
		_:
			handled = false
	if handled:
		get_viewport().set_input_as_handled()

func _spawn_at_cursor(type_name: String):
	var pid = profile_ids.get(type_name, -1)
	if pid < 0:
		return
	var pos = Vector3(randf_range(-2, 2), randf_range(-2, 2), randf_range(-2, 2))
	atom_sim.spawn_particle(pid, pos, Vector3.ZERO, Vector3.UP)

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
