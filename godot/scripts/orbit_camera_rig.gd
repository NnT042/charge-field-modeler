extends Node3D
## Orbit camera rig.
##
## Place a Camera3D as a child of this node, offset along +Z by some distance.
## This rig pivots around its own origin (the focus point). LMB drag rotates,
## MMB drag pans the focus point in the camera's local plane, and the mouse
## wheel zooms by adjusting the camera's local Z offset.

@export var min_zoom: float = 0.5
@export var max_zoom: float = 200.0
@export var zoom_step: float = 1.1
@export var rotate_sensitivity: float = 0.005
@export var pan_sensitivity: float = 0.0025
@export var pitch_limit: float = deg_to_rad(89.0)

var _camera: Camera3D
var _yaw: float = 0.0
var _pitch: float = deg_to_rad(20.0)
var _distance: float = 5.0

var _rotating: bool = false
var _panning: bool = false


func _ready() -> void:
	_camera = _find_camera()
	if _camera == null:
		push_error("OrbitCameraRig: no Camera3D child found")
		return
	_distance = _camera.position.z
	_apply_transform()


func _find_camera() -> Camera3D:
	for child in get_children():
		if child is Camera3D:
			return child
	return null


func _unhandled_input(event: InputEvent) -> void:
	if event is InputEventMouseButton:
		var mb: InputEventMouseButton = event
		match mb.button_index:
			MOUSE_BUTTON_LEFT:
				_rotating = mb.pressed
			MOUSE_BUTTON_MIDDLE:
				_panning = mb.pressed
			MOUSE_BUTTON_WHEEL_UP:
				if mb.pressed:
					_zoom(1.0 / zoom_step)
			MOUSE_BUTTON_WHEEL_DOWN:
				if mb.pressed:
					_zoom(zoom_step)
	elif event is InputEventMouseMotion:
		var mm: InputEventMouseMotion = event
		if _rotating:
			_yaw -= mm.relative.x * rotate_sensitivity
			_pitch = clamp(
				_pitch - mm.relative.y * rotate_sensitivity,
				-pitch_limit,
				pitch_limit,
			)
			_apply_transform()
		elif _panning:
			# Pan in the camera's local right/up plane, scaled by distance so
			# panning feels consistent regardless of zoom level.
			var right: Vector3 = global_transform.basis.x
			var up: Vector3 = global_transform.basis.y
			var scale: float = _distance * pan_sensitivity
			position -= right * mm.relative.x * scale
			position += up * mm.relative.y * scale


func _zoom(factor: float) -> void:
	_distance = clamp(_distance * factor, min_zoom, max_zoom)
	if _camera != null:
		_camera.position = Vector3(0.0, 0.0, _distance)


func _apply_transform() -> void:
	rotation = Vector3(_pitch, _yaw, 0.0)
	if _camera != null:
		_camera.position = Vector3(0.0, 0.0, _distance)
