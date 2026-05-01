extends Node3D
## Renders the focus particle's worldline.
##
## Three display modes cycled by T:
##   0 = OFF   — tracing disabled, nothing drawn
##   1 = LINE  — unshaded line-strip of the center point's path
##   2 = TUBE  — triangle-strip tube showing the full swept volume
##
## C clears the trace in any active mode.

@export var focus_particle_path: NodePath
@export var color: Color = Color(1.0, 0.95, 0.85, 1.0)
@export var tube_color: Color = Color(0.6, 0.55, 0.45, 0.35)
@export var tail_min_alpha: float = 0.05
@export var tube_segments: int = 10

## 0 = off, 1 = line, 2 = tube
var _display_mode: int = 1
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
	_material.cull_mode = BaseMaterial3D.CULL_DISABLED

	_mesh = ImmediateMesh.new()
	_mesh_instance = MeshInstance3D.new()
	_mesh_instance.mesh = _mesh
	add_child(_mesh_instance)

	_focus.call("set_path_enabled", _display_mode > 0)


func display_mode_label() -> String:
	match _display_mode:
		0: return "off"
		1: return "line"
		2: return "tube"
		_: return "off"


func _process(_delta: float) -> void:
	if _focus == null or _mesh == null:
		return

	_mesh.clear_surfaces()

	if _display_mode == 0:
		return

	var points: PackedVector3Array = _focus.call("get_path_points")
	if points.size() < 2:
		return

	if _display_mode == 1:
		_draw_line(points, color)
	else:
		var radii: PackedFloat32Array = _focus.call("get_path_radii")
		_draw_tube(points, radii, tube_color)


func _draw_line(points: PackedVector3Array, c: Color) -> void:
	_mesh.surface_begin(Mesh.PRIMITIVE_LINE_STRIP, _material)
	var n: int = points.size()
	var inv_last: float = 1.0 / float(n - 1)
	for i in n:
		var t: float = float(i) * inv_last
		var a: float = lerp(tail_min_alpha, 1.0, t) * c.a
		_mesh.surface_set_color(Color(c.r, c.g, c.b, a))
		_mesh.surface_add_vertex(points[i])
	_mesh.surface_end()


func _draw_tube(points: PackedVector3Array, radii: PackedFloat32Array, c: Color) -> void:
	var n: int = points.size()
	if n < 2 or radii.size() < n:
		return

	var seg: int = tube_segments
	var angle_step: float = TAU / float(seg)
	var inv_last: float = 1.0 / float(n - 1)

	# Build ring vertices for each sample point.
	# Each ring lies in the plane perpendicular to the path tangent.
	var rings: Array[PackedVector3Array] = []
	rings.resize(n)
	for i in n:
		var fwd: Vector3
		if i < n - 1:
			fwd = (points[i + 1] - points[i]).normalized()
		else:
			fwd = (points[i] - points[i - 1]).normalized()
		if fwd.is_zero_approx():
			fwd = Vector3.UP

		# Build a stable perpendicular basis.
		var up_hint: Vector3 = Vector3.UP if abs(fwd.dot(Vector3.UP)) < 0.99 else Vector3.RIGHT
		var right: Vector3 = fwd.cross(up_hint).normalized()
		var up: Vector3 = right.cross(fwd).normalized()

		var r: float = radii[i]
		var ring: PackedVector3Array = PackedVector3Array()
		ring.resize(seg)
		for j in seg:
			var theta: float = float(j) * angle_step
			ring[j] = points[i] + (right * cos(theta) + up * sin(theta)) * r
		rings[i] = ring

	_mesh.surface_begin(Mesh.PRIMITIVE_TRIANGLES, _material)
	for i in range(n - 1):
		var t0: float = float(i) * inv_last
		var t1: float = float(i + 1) * inv_last
		var a0: float = lerp(tail_min_alpha, 1.0, t0) * c.a
		var a1: float = lerp(tail_min_alpha, 1.0, t1) * c.a
		var c0: Color = Color(c.r, c.g, c.b, a0)
		var c1: Color = Color(c.r, c.g, c.b, a1)

		for j in seg:
			var j_next: int = (j + 1) % seg
			_mesh.surface_set_color(c0)
			_mesh.surface_add_vertex(rings[i][j])
			_mesh.surface_set_color(c1)
			_mesh.surface_add_vertex(rings[i + 1][j])
			_mesh.surface_set_color(c0)
			_mesh.surface_add_vertex(rings[i][j_next])

			_mesh.surface_set_color(c0)
			_mesh.surface_add_vertex(rings[i][j_next])
			_mesh.surface_set_color(c1)
			_mesh.surface_add_vertex(rings[i + 1][j])
			_mesh.surface_set_color(c1)
			_mesh.surface_add_vertex(rings[i + 1][j_next])

	_mesh.surface_end()


func cycle_display_mode() -> void:
	_display_mode = (_display_mode + 1) % 3
	if _focus != null:
		_focus.call("set_path_enabled", _display_mode > 0)
	if _display_mode == 0 and _mesh != null:
		_mesh.clear_surfaces()


func clear_trace() -> void:
	if _focus != null:
		_focus.call("clear_path")
	if _mesh != null:
		_mesh.clear_surfaces()


func _unhandled_key_input(event: InputEvent) -> void:
	if not (event is InputEventKey) or not event.pressed or event.echo:
		return
	if _focus == null:
		return
	var k: InputEventKey = event
	match k.keycode:
		KEY_T:
			cycle_display_mode()
		KEY_C:
			clear_trace()
