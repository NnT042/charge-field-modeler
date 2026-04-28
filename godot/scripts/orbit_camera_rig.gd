extends Node3D
## Orbit camera rig.
##
## Place a Camera3D as a child of this node, offset along +Z by some distance.
## This rig pivots around its own origin (the focus point). LMB drag rotates,
## MMB drag pans the focus point in the camera's local plane, and the mouse
## wheel zooms by adjusting the camera's local Z offset.
##
## Axis-snap hotkeys (regular row + numpad both work):
##   1 = Front (look down -Z, sees XY plane)
##   3 = Right (look down -X, sees YZ plane = X-spin orbit plane)
##   7 = Top   (look down -Y, sees XZ plane = Y-spin orbit plane)
##   5 = Iso   (yaw=-45°, pitch=30°)
##   0 = Reset to defaults (also resets zoom distance)

@export var min_zoom: float = 0.5
@export var max_zoom: float = 200.0
@export var zoom_step: float = 1.1
@export var rotate_sensitivity: float = 0.005
@export var pan_sensitivity: float = 0.0025
@export var pitch_limit: float = deg_to_rad(89.0)
@export var snap_duration: float = 0.25

const _DEFAULT_YAW: float = 0.0
const _DEFAULT_PITCH_DEG: float = 20.0
const _DEFAULT_DISTANCE: float = 5.0

var _camera: Camera3D
var _yaw: float = _DEFAULT_YAW
var _pitch: float = deg_to_rad(_DEFAULT_PITCH_DEG)
var _distance: float = _DEFAULT_DISTANCE

var _rotating: bool = false
var _panning: bool = false
var _snap_tween: Tween


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
				if mb.pressed:
					_kill_snap()
			MOUSE_BUTTON_MIDDLE:
				_panning = mb.pressed
			MOUSE_BUTTON_WHEEL_UP:
				if mb.pressed:
					_kill_snap()
					_zoom(1.0 / zoom_step)
			MOUSE_BUTTON_WHEEL_DOWN:
				if mb.pressed:
					_kill_snap()
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


func _unhandled_key_input(event: InputEvent) -> void:
	if not (event is InputEventKey):
		return
	var k: InputEventKey = event
	if not k.pressed or k.echo:
		return
	var handled: bool = true
	match k.keycode:
		KEY_1, KEY_KP_1:
			_snap_to(0.0, 0.0)
		KEY_3, KEY_KP_3:
			_snap_to(-PI * 0.5, 0.0)
		KEY_7, KEY_KP_7:
			_snap_to(0.0, pitch_limit)
		KEY_5, KEY_KP_5:
			_snap_to(-PI * 0.25, deg_to_rad(30.0))
		KEY_0, KEY_KP_0:
			_snap_to(_DEFAULT_YAW, deg_to_rad(_DEFAULT_PITCH_DEG), _DEFAULT_DISTANCE)
		_:
			handled = false
	if handled:
		get_viewport().set_input_as_handled()


func _zoom(factor: float) -> void:
	_distance = clamp(_distance * factor, min_zoom, max_zoom)
	if _camera != null:
		_camera.position = Vector3(0.0, 0.0, _distance)


func _apply_transform() -> void:
	rotation = Vector3(_pitch, _yaw, 0.0)
	if _camera != null:
		_camera.position = Vector3(0.0, 0.0, _distance)


# Tween yaw/pitch (and optionally distance) to a target view. -1.0 distance
# means "leave zoom alone". Yaw uses shortest-path wrap so a snap from yaw=0
# to yaw=-π/2 doesn't go the long way around.
func _snap_to(target_yaw: float, target_pitch: float, target_distance: float = -1.0) -> void:
	_kill_snap()
	var resolved_yaw: float = _yaw + wrapf(target_yaw - _yaw, -PI, PI)
	var resolved_distance: float = target_distance if target_distance > 0.0 else _distance
	var from_state: Vector3 = Vector3(_yaw, _pitch, _distance)
	var to_state: Vector3 = Vector3(resolved_yaw, target_pitch, resolved_distance)
	_snap_tween = create_tween().set_trans(Tween.TRANS_SINE).set_ease(Tween.EASE_IN_OUT)
	_snap_tween.tween_method(_set_orbit_state, from_state, to_state, snap_duration)


func _set_orbit_state(state: Vector3) -> void:
	_yaw = state.x
	_pitch = clamp(state.y, -pitch_limit, pitch_limit)
	_distance = clamp(state.z, min_zoom, max_zoom)
	_apply_transform()


func _kill_snap() -> void:
	if _snap_tween != null and _snap_tween.is_valid():
		_snap_tween.kill()
	_snap_tween = null
