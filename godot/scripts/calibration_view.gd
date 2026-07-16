extends Node3D
## Eyeball renderer for the Field Calibration Mode baked-loop hitbox (CM-0 dev
## scene). Draws the CalibrationView node's live loop as a closed line strip
## plus a pole-axis segment, and forwards keyboard controls to the Rust node.
##
## Controls:
##   Z / X / C        proton / neutron / electron
##   Up / Down        outer_spin (swing rate) +/-
##   Left / Right     tumble rate +/-
##   Space            nudge (one scaled random impulse)
##   G                ambient field on / off (watch outer_spin pump)
##   T                cycle temperature preset (cold_void→room→hot→solar_core)
##   B                field balanced (0.5) vs photon-rich (0.67, spins up)
##   V                field directional (+X drift) vs isotropic
##   P                pause / resume
##   R                reset current particle
##   (camera: mouse-drag orbit, wheel zoom, F fit, 5 iso, 0 default)

@export var view_path: NodePath
@export var label_path: NodePath
@export var loop_color: Color = Color(0.55, 1.0, 0.95, 0.95)
@export var pole_color: Color = Color(1.0, 0.35, 0.2, 0.9)

var _view: Node3D
var _label: Label
var _mesh: ArrayMesh
var _mesh_instance: MeshInstance3D
var _material: StandardMaterial3D


func _ready() -> void:
	_view = get_node_or_null(view_path) as Node3D
	if _view == null:
		push_error("calibration_view.gd: view_path does not resolve to the CalibrationView node")
		return
	_label = get_node_or_null(label_path) as Label

	_material = StandardMaterial3D.new()
	_material.shading_mode = StandardMaterial3D.SHADING_MODE_UNSHADED
	_material.vertex_color_use_as_albedo = true
	_material.transparency = StandardMaterial3D.TRANSPARENCY_ALPHA
	_material.disable_receive_shadows = true
	_material.cull_mode = BaseMaterial3D.CULL_DISABLED

	_mesh = ArrayMesh.new()
	_mesh_instance = MeshInstance3D.new()
	_mesh_instance.mesh = _mesh
	add_child(_mesh_instance)


func _process(_delta: float) -> void:
	if _view == null or _mesh == null:
		return
	_mesh.clear_surfaces()

	var loop: PackedVector3Array = _view.call("get_loop_points")
	if loop.size() >= 2:
		# Close the loop so the swept curve reads as a ring, not an open strip.
		loop.push_back(loop[0])
		_add_surface(Mesh.PRIMITIVE_LINE_STRIP, loop, _solid_colors(loop.size(), loop_color))

	var pole: PackedVector3Array = _view.call("get_pole_segment")
	if pole.size() == 2:
		_add_surface(Mesh.PRIMITIVE_LINES, pole, _solid_colors(2, pole_color))

	if _label != null:
		_label.text = str(_view.call("get_status_text"))


func _solid_colors(n: int, c: Color) -> PackedColorArray:
	var colors: PackedColorArray = PackedColorArray()
	colors.resize(n)
	for i in n:
		colors[i] = c
	return colors


func _add_surface(primitive: Mesh.PrimitiveType, verts: PackedVector3Array, colors: PackedColorArray) -> void:
	var arrays: Array = []
	arrays.resize(Mesh.ARRAY_MAX)
	arrays[Mesh.ARRAY_VERTEX] = verts
	arrays[Mesh.ARRAY_COLOR] = colors
	_mesh.add_surface_from_arrays(primitive, arrays)
	_mesh.surface_set_material(_mesh.get_surface_count() - 1, _material)


func _unhandled_key_input(event: InputEvent) -> void:
	if not (event is InputEventKey) or _view == null:
		return
	var k: InputEventKey = event
	if not k.pressed or k.echo:
		return
	var handled: bool = true
	match k.keycode:
		KEY_Z:
			_view.call("set_kind", 0)
		KEY_X:
			_view.call("set_kind", 1)
		KEY_C:
			_view.call("set_kind", 2)
		KEY_UP:
			_view.call("bump_outer_spin", 0.05)
		KEY_DOWN:
			_view.call("bump_outer_spin", -0.05)
		KEY_LEFT:
			_view.call("bump_tumble", -0.05)
		KEY_RIGHT:
			_view.call("bump_tumble", 0.05)
		KEY_SPACE:
			_view.call("nudge")
		KEY_P:
			_view.call("toggle_pause")
		KEY_R:
			_view.call("reset")
		KEY_G:
			_view.call("toggle_field")
		KEY_T:
			_view.call("cycle_preset")
		KEY_B:
			_view.call("toggle_balance")
		KEY_V:
			_view.call("toggle_directional")
		_:
			handled = false
	if handled:
		get_viewport().set_input_as_handled()
