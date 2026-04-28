extends Node3D
## Transparent orbital-shell reference for the focus particle.
##
## Center: tracks the outermost level's orbit center (the inner levels'
## accumulated position) so the particle always lies on the ghost sphere's
## equatorial ring. This matches Wheeler's wireframe-ball visualization.
##
## Radius: effective_radius × WORLD_SCALE (0.5). At axial level the ghost
## sphere is the same size as the photon mesh (both radius 0.5 world units).
##
## Orientation: the ghost sphere's +Y pole aligns with the outermost spin's
## rotation axis, so the equatorial ring marks the orbital plane and the
## poles mark the rotation axis — readable at a glance without labels.
##
## Hotkey: G toggles visibility.

const WORLD_SCALE: float = 0.5

@export var focus_particle_path: NodePath
@export var visible_by_default: bool = true

var _focus: Node3D
var _mesh_instance: MeshInstance3D


func _ready() -> void:
	_focus = get_node_or_null(focus_particle_path) as Node3D
	if _focus == null:
		push_error("GhostSphere: focus_particle_path does not resolve to a Node3D")
		return

	var mat := ShaderMaterial.new()
	mat.shader = load("res://shaders/visual/ghost_sphere.gdshader")

	var mesh := SphereMesh.new()
	mesh.radius = 1.0
	mesh.height = 2.0
	mesh.radial_segments = 64
	mesh.rings = 32

	_mesh_instance = MeshInstance3D.new()
	_mesh_instance.mesh = mesh
	_mesh_instance.material_override = mat
	_mesh_instance.cast_shadow = GeometryInstance3D.SHADOW_CASTING_SETTING_OFF
	add_child(_mesh_instance)

	visible = visible_by_default


func _process(_delta: float) -> void:
	if _focus == null or _mesh_instance == null:
		return

	# Ghost sphere is fixed at world origin — it marks the outermost level's
	# orbital shell, not the particle's current position. Only scale and
	# orientation change.
	var r: float = float(_focus.call("effective_radius")) * WORLD_SCALE
	scale = Vector3(r, r, r)

	# Orient so +Y pole aligns with the outermost spin's rotation axis.
	var count: int = int(_focus.call("level_count"))
	var label: String = String(_focus.call("spin_type_label", count))
	_align_pole_to(_label_to_axis(label))


func _unhandled_key_input(event: InputEvent) -> void:
	if not (event is InputEventKey) or not event.pressed or event.echo:
		return
	if (event as InputEventKey).keycode == KEY_G:
		visible = not visible


# Rotate this node so its local +Y axis points to `axis` in world space.
func _align_pole_to(axis: Vector3) -> void:
	var up := Vector3.UP
	if axis.is_equal_approx(up):
		quaternion = Quaternion.IDENTITY
	elif axis.is_equal_approx(-up):
		quaternion = Quaternion(Vector3.RIGHT, PI)
	else:
		quaternion = Quaternion(up.cross(axis).normalized(), up.angle_to(axis))


# Map the Rust spin-type label to the rotation axis in world space.
func _label_to_axis(label: String) -> Vector3:
	match label:
		"x": return Vector3.RIGHT
		"z": return Vector3(0.0, 0.0, 1.0)
		_:   return Vector3.UP  # "axial" and "y" both rotate around world Y
