extends Node3D
## Renders the focus particle's worldline as an unshaded line-strip with a
## per-vertex alpha fade (oldest = transparent, newest = opaque). Reads the
## points from the FocusParticle Rust class each frame and rebuilds an
## ImmediateMesh.
##
## Hotkeys: T toggles tracing, C clears the trace.

@export var focus_particle_path: NodePath
@export var color: Color = Color(1.0, 0.95, 0.85, 1.0)
@export var tail_min_alpha: float = 0.05

var _focus: Node3D
var _mesh_instance: MeshInstance3D
var _mesh: ImmediateMesh
var _material: StandardMaterial3D


func _ready() -> void:
	_focus = get_node_or_null(focus_particle_path) as Node3D
	if _focus == null:
		push_error("PathTraceLine: focus_particle_path does not resolve to a Node3D")
		return

	_material = StandardMaterial3D.new()
	_material.shading_mode = StandardMaterial3D.SHADING_MODE_UNSHADED
	_material.vertex_color_use_as_albedo = true
	_material.transparency = StandardMaterial3D.TRANSPARENCY_ALPHA
	_material.disable_receive_shadows = true

	_mesh = ImmediateMesh.new()
	_mesh_instance = MeshInstance3D.new()
	_mesh_instance.mesh = _mesh
	add_child(_mesh_instance)


func _process(_delta: float) -> void:
	if _focus == null or _mesh == null:
		return

	var points: PackedVector3Array = _focus.call("get_path_points")
	_mesh.clear_surfaces()
	if points.size() < 2:
		return

	_mesh.surface_begin(Mesh.PRIMITIVE_LINE_STRIP, _material)
	var n: int = points.size()
	var inv_last: float = 1.0 / float(n - 1)
	for i in n:
		# i=0 is oldest (faded), i=n-1 is newest (fully opaque).
		var t: float = float(i) * inv_last
		var a: float = lerp(tail_min_alpha, 1.0, t) * color.a
		_mesh.surface_set_color(Color(color.r, color.g, color.b, a))
		_mesh.surface_add_vertex(points[i])
	_mesh.surface_end()


func _unhandled_key_input(event: InputEvent) -> void:
	if not (event is InputEventKey) or not event.pressed or event.echo:
		return
	if _focus == null:
		return
	var k: InputEventKey = event
	match k.keycode:
		KEY_T:
			var enabled: bool = bool(_focus.call("is_path_enabled"))
			_focus.call("set_path_enabled", not enabled)
			if enabled:
				# Was on, now off — drop the visible strip immediately.
				if _mesh != null:
					_mesh.clear_surfaces()
		KEY_C:
			_focus.call("clear_path")
			if _mesh != null:
				_mesh.clear_surfaces()
