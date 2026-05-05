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
	# Tier 1 — Photon: classic RGB
	"x1": Color(1.0, 0.35, 0.35),
	"y1": Color(0.35, 1.0, 0.35),
	"z1": Color(0.4, 0.6, 1.0),
	# Tier 2 — Electron: warm shifted (axial = white-blue to mark new precession axis)
	"axial2": Color(0.7, 0.8, 1.0),
	"x2": Color(1.0, 0.6, 0.15),
	"y2": Color(0.15, 0.85, 0.85),
	"z2": Color(0.85, 0.3, 0.85),
	# Tier 3 — Baryon: earth tones (axial = pale gold)
	"axial3": Color(0.9, 0.85, 0.6),
	"x3": Color(0.95, 0.8, 0.2),
	"y3": Color(0.2, 0.7, 0.6),
	"z3": Color(0.55, 0.3, 0.9),
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
		var orient: Quaternion = Quaternion(_focus.call("ghost_sphere_orientation", i))

		mi.position = center
		mi.scale = Vector3(r, r, r)
		mi.material_override = _materials.get(label, _materials["x1"])
		var base_axis: Vector3 = _label_to_axis(label)
		var base_quat: Quaternion = _axis_to_quat(base_axis)
		mi.quaternion = orient * base_quat


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


func _axis_to_quat(axis: Vector3) -> Quaternion:
	var up := Vector3.UP
	if axis.is_equal_approx(up):
		return Quaternion.IDENTITY
	elif axis.is_equal_approx(-up):
		return Quaternion(Vector3.RIGHT, PI)
	else:
		return Quaternion(up.cross(axis).normalized(), up.angle_to(axis))


func _label_to_axis(label: String) -> Vector3:
	if label.begins_with("axial"):
		return Vector3(0.0, 0.0, 1.0)
	elif label.begins_with("x"):
		return Vector3.RIGHT
	elif label.begins_with("z"):
		return Vector3(0.0, 0.0, 1.0)
	else:
		return Vector3.UP
