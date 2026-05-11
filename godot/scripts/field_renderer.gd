extends Node3D

@export var focus_particle_path: NodePath
@export var default_photon_count: int = 500000

var _quad_mesh: QuadMesh
var _multimesh: MultiMesh
var _multimesh_instance: MultiMeshInstance3D


func _ready() -> void:
	var shader: Shader = load("res://shaders/visual/field_dot.gdshader")
	var material := ShaderMaterial.new()
	material.shader = shader

	_quad_mesh = QuadMesh.new()
	_quad_mesh.size = Vector2(0.5, 0.5)
	_quad_mesh.material = material

	_multimesh = MultiMesh.new()
	_multimesh.transform_format = MultiMesh.TRANSFORM_3D
	_multimesh.use_colors = true
	_multimesh.mesh = _quad_mesh
	_multimesh.instance_count = 0

	_multimesh_instance = MultiMeshInstance3D.new()
	_multimesh_instance.multimesh = _multimesh
	add_child(_multimesh_instance)


func set_instance_count(count: int) -> void:
	_multimesh.instance_count = count


func update_from_buffer(buffer: PackedFloat32Array) -> void:
	_multimesh.set_buffer(buffer)


func set_visible_field(vis: bool) -> void:
	_multimesh_instance.visible = vis


func toggle_visible() -> bool:
	var new_state: bool = !_multimesh_instance.visible
	_multimesh_instance.visible = new_state
	return new_state


func get_photon_count() -> int:
	return _multimesh.instance_count


func get_multimesh() -> MultiMesh:
	return _multimesh