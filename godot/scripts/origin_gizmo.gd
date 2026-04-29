extends Node3D
## Static visual reference at world origin showing the three coordinate axes
## as colored unshaded line segments (X=red, Y=green, Z=blue) with tick marks
## at each natural-unit interval (1r = 0.5 world units).

const WORLD_SCALE: float = 0.5
const TICK_SIZE: float = 0.04
const AXIS_LENGTH: float = 4.5


func _ready() -> void:
	var mesh: ImmediateMesh = ImmediateMesh.new()
	var material: StandardMaterial3D = StandardMaterial3D.new()
	material.shading_mode = StandardMaterial3D.SHADING_MODE_UNSHADED
	material.vertex_color_use_as_albedo = true
	material.disable_receive_shadows = true

	var x_color := Color(1.0, 0.35, 0.35)
	var y_color := Color(0.35, 1.0, 0.35)
	var z_color := Color(0.4, 0.6, 1.0)
	var tick_color := Color(0.7, 0.7, 0.7, 0.6)

	mesh.surface_begin(Mesh.PRIMITIVE_LINES, material)

	# +X axis
	mesh.surface_set_color(x_color)
	mesh.surface_add_vertex(Vector3.ZERO)
	mesh.surface_add_vertex(Vector3(AXIS_LENGTH, 0.0, 0.0))

	# +Y axis
	mesh.surface_set_color(y_color)
	mesh.surface_add_vertex(Vector3.ZERO)
	mesh.surface_add_vertex(Vector3(0.0, AXIS_LENGTH, 0.0))

	# +Z axis
	mesh.surface_set_color(z_color)
	mesh.surface_add_vertex(Vector3.ZERO)
	mesh.surface_add_vertex(Vector3(0.0, 0.0, AXIS_LENGTH))

	# Tick marks at each natural-unit interval along each axis.
	# Key amplitudes: 1r, 2r, 4r, 8r (levels 1-4).
	var r: int = 1
	while float(r) * WORLD_SCALE < AXIS_LENGTH:
		var d: float = float(r) * WORLD_SCALE
		var big: bool = (r == 1 or r == 2 or r == 4 or r == 8)
		var t: float = TICK_SIZE * (2.0 if big else 1.0)

		# X-axis ticks (perpendicular in Y)
		mesh.surface_set_color(x_color if big else tick_color)
		mesh.surface_add_vertex(Vector3(d, -t, 0.0))
		mesh.surface_add_vertex(Vector3(d, t, 0.0))

		# Y-axis ticks (perpendicular in X)
		mesh.surface_set_color(y_color if big else tick_color)
		mesh.surface_add_vertex(Vector3(-t, d, 0.0))
		mesh.surface_add_vertex(Vector3(t, d, 0.0))

		# Z-axis ticks (perpendicular in Y)
		mesh.surface_set_color(z_color if big else tick_color)
		mesh.surface_add_vertex(Vector3(0.0, -t, d))
		mesh.surface_add_vertex(Vector3(0.0, t, d))

		r += 1

	mesh.surface_end()

	var mesh_instance: MeshInstance3D = MeshInstance3D.new()
	mesh_instance.mesh = mesh
	add_child(mesh_instance)
