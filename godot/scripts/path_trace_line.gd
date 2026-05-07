extends Node3D
## Renders the focus particle's worldline.
##
## Five display modes cycled by T:
##   0 = OFF     — tracing disabled, nothing drawn
##   1 = LINE    — unshaded line-strip of the center point's path
##   2 = TUBE    — triangle-strip tube showing the full swept volume
##   3 = SURFACE — line-strip following a point on the particle's pole
##   4 = SNAKE   — fading tube tail showing 1/4 rotation of the outermost ghost
##
## C clears the trace in any active mode.

@export var focus_particle_path: NodePath
@export var color: Color = Color(1.0, 0.95, 0.85, 1.0)
@export var tube_color: Color = Color(0.6, 0.55, 0.45, 0.35)
@export var surface_color: Color = Color(0.4, 0.8, 1.0, 0.9)
@export var crest_color: Color = Color(1.0, 0.4, 0.1, 0.9)
@export var trough_color: Color = Color(0.2, 0.6, 1.0, 0.9)
@export var marker_size: float = 0.1
@export var tail_min_alpha: float = 0.05
@export var tube_segments: int = 10

## 0 = off, 1 = line, 2 = tube, 3 = surface, 4 = snake
var _display_mode: int = 1
var _show_marker_lines: bool = false
var _focus: Node3D
var _mesh_instance: MeshInstance3D
var _mesh: ArrayMesh
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

	_mesh = ArrayMesh.new()
	_mesh_instance = MeshInstance3D.new()
	_mesh_instance.mesh = _mesh
	add_child(_mesh_instance)

	_focus.call("set_path_enabled", _display_mode > 0)


func display_mode_label() -> String:
	match _display_mode:
		0: return "off"
		1: return "line"
		2: return "tube"
		3: return "surface"
		4: return "snake"
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

	var spectral: Color = Color(_focus.call("get_wavelength_color"))
	var in_visible: bool = not spectral.is_equal_approx(Color(0.5, 0.5, 0.5))
	var line_c: Color = Color(spectral.r, spectral.g, spectral.b, color.a) if in_visible else color
	var tube_c: Color = Color(spectral.r, spectral.g, spectral.b, tube_color.a) if in_visible else tube_color
	var surf_c: Color = Color(spectral.r, spectral.g, spectral.b, surface_color.a) if in_visible else surface_color

	if _display_mode == 1:
		_draw_line(points, line_c)
	elif _display_mode == 2:
		_draw_tube_from_rust(tube_c, -1)
	elif _display_mode == 3:
		var surface_pts: PackedVector3Array = _focus.call("get_surface_path_points")
		if surface_pts.size() >= 2:
			_draw_line(surface_pts, surf_c)
	elif _display_mode == 4:
		var snake_len: int = int(_focus.call("snake_trace_length"))
		_draw_tube_from_rust(tube_c, snake_len)

	if _display_mode > 0:
		var crests: PackedVector3Array = _focus.call("get_crest_markers")
		var troughs: PackedVector3Array = _focus.call("get_trough_markers")
		if crests.size() > 0:
			_draw_markers(crests, crest_color)
		if troughs.size() > 0:
			_draw_markers(troughs, trough_color)
		if _show_marker_lines:
			if crests.size() > 1:
				_draw_marker_chain(crests, Color(crest_color.r, crest_color.g, crest_color.b, 0.5))
			if troughs.size() > 1:
				_draw_marker_chain(troughs, Color(trough_color.r, trough_color.g, trough_color.b, 0.5))


func _add_surface(primitive: Mesh.PrimitiveType, verts: PackedVector3Array, colors: PackedColorArray) -> void:
	var arrays: Array = []
	arrays.resize(Mesh.ARRAY_MAX)
	arrays[Mesh.ARRAY_VERTEX] = verts
	arrays[Mesh.ARRAY_COLOR] = colors
	_mesh.add_surface_from_arrays(primitive, arrays)
	_mesh.surface_set_material(_mesh.get_surface_count() - 1, _material)


func _draw_line(points: PackedVector3Array, c: Color) -> void:
	var n: int = points.size()
	var colors: PackedColorArray = PackedColorArray()
	colors.resize(n)
	var inv_last: float = 1.0 / float(n - 1)
	for i in n:
		var t: float = float(i) * inv_last
		var a: float = lerp(tail_min_alpha, 1.0, t) * c.a
		colors[i] = Color(c.r, c.g, c.b, a)
	_add_surface(Mesh.PRIMITIVE_LINE_STRIP, points, colors)


func _draw_tube_from_rust(c: Color, max_samples: int) -> void:
	var count: int = int(_focus.call("build_tube_mesh", tube_segments, max_samples,
		c.r, c.g, c.b, c.a, tail_min_alpha))
	if count > 0:
		_add_surface(Mesh.PRIMITIVE_TRIANGLES,
			_focus.call("get_tube_vertices"),
			_focus.call("get_tube_colors"))


func _draw_markers(points: PackedVector3Array, c: Color) -> void:
	var n: int = points.size()
	var s: float = marker_size
	var verts: PackedVector3Array = PackedVector3Array()
	verts.resize(n * 6)
	var colors: PackedColorArray = PackedColorArray()
	colors.resize(n * 6)
	for i in n:
		var p: Vector3 = points[i]
		var base: int = i * 6
		verts[base] = p + Vector3(s, 0, 0); colors[base] = c
		verts[base + 1] = p - Vector3(s, 0, 0); colors[base + 1] = c
		verts[base + 2] = p + Vector3(0, s, 0); colors[base + 2] = c
		verts[base + 3] = p - Vector3(0, s, 0); colors[base + 3] = c
		verts[base + 4] = p + Vector3(0, 0, s); colors[base + 4] = c
		verts[base + 5] = p - Vector3(0, 0, s); colors[base + 5] = c
	_add_surface(Mesh.PRIMITIVE_LINES, verts, colors)


func _draw_marker_chain(points: PackedVector3Array, c: Color) -> void:
	if points.size() < 2:
		return
	var pair_count: int = points.size() / 2
	var verts: PackedVector3Array = PackedVector3Array()
	verts.resize(pair_count * 2)
	var colors: PackedColorArray = PackedColorArray()
	colors.resize(pair_count * 2)
	for i in pair_count:
		var base: int = i * 2
		verts[base] = points[base]; colors[base] = c
		verts[base + 1] = points[base + 1]; colors[base + 1] = c
	_add_surface(Mesh.PRIMITIVE_LINES, verts, colors)


func cycle_display_mode() -> void:
	_display_mode = (_display_mode + 1) % 5
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
		KEY_M:
			_show_marker_lines = not _show_marker_lines
