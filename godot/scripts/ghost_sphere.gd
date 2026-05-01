extends Node3D
## Renders one transparent ghost sphere per active orbital level.
##
## In the hierarchical nesting model the outermost ghost sits at origin;
## inner ghosts orbit inside it, carried by the outer levels' rotations.
## Each ghost is color-coded to its spin axis: X=red, Y=green, Z=blue.
##
## Hotkey: G cycles display mode: off → all → outermost only → off.

@export var focus_particle_path: NodePath
@export var visible_by_default: bool = true

## 0 = off, 1 = all levels, 2 = outermost only
var _display_mode: int = 1
var _focus: Node3D
var _meshes: Array[MeshInstance3D] = []
var _shader: Shader
var _materials: Dictionary = {}

const AXIS_COLORS: Dictionary = {
	"x": Color(1.0, 0.35, 0.35),
	"y": Color(0.35, 1.0, 0.35),
	"z": Color(0.4, 0.6, 1.0),
}


func _ready() -> void:
	_focus = get_node_or_null(focus_particle_path) as Node3D
	if _focus == null:
		push_error("GhostSphere: focus_particle_path does not resolve to a Node3D")
		return

	_shader = load("res://shaders/visual/ghost_sphere.gdshader")

	for label in AXIS_COLORS:
		var mat := ShaderMaterial.new()
		mat.shader = _shader
		mat.set_shader_parameter("ghost_color", AXIS_COLORS[label])
		_materials[label] = mat

	_display_mode = 1 if visible_by_default else 0
	visible = _display_mode > 0


func _process(_delta: float) -> void:
	if _focus == null or _display_mode == 0:
		_sync_mesh_count(0)
		return

	var total: int = int(_focus.call("ghost_sphere_count"))
	var first: int = total - 1 if _display_mode == 2 and total > 0 else 0
	var show_count: int = total - first
	_sync_mesh_count(show_count)

	for j in range(show_count):
		var i: int = first + j
		var mi: MeshInstance3D = _meshes[j]
		var center: Vector3 = Vector3(_focus.call("ghost_sphere_center", i))
		var r: float = float(_focus.call("ghost_sphere_radius", i))
		var label: String = String(_focus.call("ghost_sphere_spin_label", i))

		mi.position = center
		mi.scale = Vector3(r, r, r)
		mi.material_override = _materials.get(label, _materials["x"])
		_align_mesh_to(mi, _label_to_axis(label))


func _sync_mesh_count(target: int) -> void:
	while _meshes.size() < target:
		var mesh := SphereMesh.new()
		mesh.radius = 1.0
		mesh.height = 2.0
		mesh.radial_segments = 64
		mesh.rings = 32

		var mi := MeshInstance3D.new()
		mi.mesh = mesh
		mi.cast_shadow = GeometryInstance3D.SHADOW_CASTING_SETTING_OFF
		add_child(mi)
		_meshes.append(mi)

	while _meshes.size() > target:
		var mi: MeshInstance3D = _meshes.pop_back()
		mi.queue_free()


func display_mode_label() -> String:
	match _display_mode:
		0: return "off"
		1: return "all"
		2: return "outer"
		_: return "off"


func cycle_display_mode() -> void:
	_display_mode = (_display_mode + 1) % 3
	visible = _display_mode > 0


func _unhandled_key_input(event: InputEvent) -> void:
	if not (event is InputEventKey) or not event.pressed or event.echo:
		return
	if (event as InputEventKey).keycode == KEY_G:
		cycle_display_mode()


func _align_mesh_to(mi: MeshInstance3D, axis: Vector3) -> void:
	var up := Vector3.UP
	if axis.is_equal_approx(up):
		mi.quaternion = Quaternion.IDENTITY
	elif axis.is_equal_approx(-up):
		mi.quaternion = Quaternion(Vector3.RIGHT, PI)
	else:
		mi.quaternion = Quaternion(up.cross(axis).normalized(), up.angle_to(axis))


func _label_to_axis(label: String) -> Vector3:
	match label:
		"x": return Vector3.RIGHT
		"z": return Vector3(0.0, 0.0, 1.0)
		_:   return Vector3.UP
