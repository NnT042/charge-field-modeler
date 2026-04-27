extends CanvasLayer
## M2 HUD: drives the spin stack (levels 1-4, photon tier).
##
## - Each level row has a slider (−c..+c) plus an activate button. The button
##   is visible only for the next-to-activate level, and only once the
##   currently-topmost level has saturated to ±c.
## - The time-scale slider is logarithmic (10^x). 0 = unity → TAU rad/s at c
##   (one revolution per second, matching M1's feel). Negative = slow-mo.
## - Readouts mirror the Rust side (signature, classification, effective
##   radius, top-level ω, FPS) and update each frame.

const PHASE_MAX_LEVEL: int = 4  ## M2 covers the photon tier (levels 1..4).
const BASE_TIME_SCALE: float = TAU  ## Time scale at slider position 0.0.

@export var focus_particle_path: NodePath

var _focus: Node3D
var _level_labels: Array[Label] = []
var _sliders: Array[HSlider] = []
var _value_labels: Array[Label] = []
var _activate_buttons: Array[Button] = []


func _ready() -> void:
	_focus = get_node_or_null(focus_particle_path) as Node3D
	if _focus == null:
		push_error("HUD: focus_particle_path does not resolve to a Node3D")
		return

	_level_labels = [%Lvl1Label, %Lvl2Label, %Lvl3Label, %Lvl4Label]
	_sliders = [%Lvl1Slider, %Lvl2Slider, %Lvl3Slider, %Lvl4Slider]
	_value_labels = [%Lvl1Value, %Lvl2Value, %Lvl3Value, %Lvl4Value]
	_activate_buttons = [%Lvl1Activate, %Lvl2Activate, %Lvl3Activate, %Lvl4Activate]

	for i in PHASE_MAX_LEVEL:
		var level: int = i + 1
		_level_labels[i].text = _format_level_label(level)
		_sliders[i].value_changed.connect(_on_slider_changed.bind(level))
		_activate_buttons[i].pressed.connect(_on_activate_pressed.bind(level))

	%TimeScaleSlider.value_changed.connect(_on_time_scale_changed)
	_apply_time_scale(%TimeScaleSlider.value)
	_refresh_row_states()


func _process(_delta: float) -> void:
	if _focus == null:
		return

	var active: int = int(_focus.call("level_count"))
	for i in PHASE_MAX_LEVEL:
		var level: int = i + 1
		if level <= active:
			var v: float = float(_focus.call("get_level_velocity", level))
			if not _sliders[i].has_focus():
				_sliders[i].set_value_no_signal(v)
			_value_labels[i].text = "%+0.3f c" % v
		else:
			_value_labels[i].text = "—"

	_refresh_row_states()

	%ParticleLabel.text = String(_focus.call("classification"))
	%SignatureLabel.text = String(_focus.call("signature"))
	%RadiusLabel.text = "%d r" % int(_focus.call("effective_radius"))
	if active > 0:
		var top_v: float = float(_focus.call("get_level_velocity", active))
		%TopOmegaLabel.text = "%+0.3f c (lvl %d)" % [top_v, active]
	else:
		%TopOmegaLabel.text = "—"
	%FpsLabel.text = "%d" % Engine.get_frames_per_second()


func _on_slider_changed(value: float, level: int) -> void:
	if _focus == null:
		return
	_focus.call("set_level_velocity", level, value)


func _on_activate_pressed(level: int) -> void:
	if _focus == null:
		return
	# Walk activate_next() up to the requested level (defensive — in practice
	# the button is only enabled for level == current_top + 1).
	while int(_focus.call("level_count")) < level:
		var added: int = int(_focus.call("activate_next"))
		if added == 0:
			break
	_refresh_row_states()


func _on_time_scale_changed(slider_log: float) -> void:
	_apply_time_scale(slider_log)


func _apply_time_scale(slider_log: float) -> void:
	var factor: float = pow(10.0, slider_log)
	var ts: float = BASE_TIME_SCALE * factor
	if _focus != null:
		_focus.call("set_time_scale", ts)
	%TimeScaleValue.text = "%.3f× (%.2f rad/s @ c)" % [factor, ts]


func _refresh_row_states() -> void:
	if _focus == null:
		return
	var active: int = int(_focus.call("level_count"))
	var can_next: bool = bool(_focus.call("can_activate_next"))
	for i in PHASE_MAX_LEVEL:
		var level: int = i + 1
		var is_active: bool = level <= active
		_sliders[i].editable = is_active
		_sliders[i].modulate = Color(1, 1, 1, 1) if is_active else Color(1, 1, 1, 0.35)
		var is_next: bool = (level == active + 1) and can_next
		_activate_buttons[i].visible = is_next


func _format_level_label(level: int) -> String:
	if _focus == null:
		return "Lvl %d" % level
	var spin_type: String = String(_focus.call("spin_type_label", level))
	var tier: String = String(_focus.call("tier_label", level))
	var amp: float = float(_focus.call("level_amplitude", level))
	return "Lvl %d — %s (%s, r=%d)" % [level, spin_type, tier, int(amp)]
