extends Node3D
## Static visual reference at world origin showing the three coordinate axes
## as colored unshaded line segments (X=red, Y=green, Z=blue).


func _ready() -> void:
	var mesh: ImmediateMesh = ImmediateMesh.new()
	var material: StandardMaterial3D = StandardMaterial3D.new()
	material.shading_mode = StandardMaterial3D.SHADING_MODE_UNSHADED
	material.vertex_color_use_as_albedo = true
	material.disable_receive_shadows = true

	mesh.surface_begin(Mesh.PRIMITIVE_LINES, material)

	# +X axis
	mesh.surface_set_color(Color(1.0, 0.35, 0.35))
	mesh.surface_add_vertex(Vector3.ZERO)
	mesh.surface_add_vertex(Vector3(2.5, 0.0, 0.0))

	# +Y axis
	mesh.surface_set_color(Color(0.35, 1.0, 0.35))
	mesh.surface_add_vertex(Vector3.ZERO)
	mesh.surface_add_vertex(Vector3(0.0, 2.5, 0.0))

	# +Z axis
	mesh.surface_set_color(Color(0.4, 0.6, 1.0))
	mesh.surface_add_vertex(Vector3.ZERO)
	mesh.surface_add_vertex(Vector3(0.0, 0.0, 2.5))

	mesh.surface_end()

	var mesh_instance: MeshInstance3D = MeshInstance3D.new()
	mesh_instance.mesh = mesh
	add_child(mesh_instance)
