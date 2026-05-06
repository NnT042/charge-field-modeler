extends CanvasLayer
## M4 HUD: toolbar + tier-tabbed slider panel for the spin stack (levels 1-16).

const PHASE_MAX_LEVEL: int = 16
const LEVELS_PER_TIER: int = 4
const TIER_COUNT: int = 4
const BASE_TIME_SCALE: float = TAU

@export var focus_particle_path: NodePath
@export var camera_rig_path: NodePath
@export var ghost_sphere_path: NodePath
@export var path_trace_path: NodePath

var _focus: Node3D
var _camera_rig: Node3D
var _ghost: Node3D
var _trace: Node3D
var _tab_container: TabContainer
var _level_labels: Array[Label] = []
var _sliders: Array[HSlider] = []
var _value_labels: Array[Label] = []
var _reverse_buttons: Array[Button] = []
var _paused: bool = false
var _prev_level_count: int = 0


func _ready() -> void:
	_focus = get_node_or_null(focus_particle_path) as Node3D
	_camera_rig = get_node_or_null(camera_rig_path) as Node3D
	_ghost = get_node_or_null(ghost_sphere_path) as Node3D
	_trace = get_node_or_null(path_trace_path) as Node3D

	if _focus == null:
		push_error("HUD: focus_particle_path does not resolve to a Node3D")
		return

	_build_tier_tabs()

	for i in PHASE_MAX_LEVEL:
		var level: int = i + 1
		_sliders[i].value_changed.connect(_on_slider_changed.bind(level))
		_reverse_buttons[i].pressed.connect(_on_reverse_pressed.bind(level))

	%TimeScaleSlider.value_changed.connect(_on_time_scale_changed)
	%LinearSpeedSlider.value_changed.connect(_on_linear_speed_changed)
	_apply_time_scale(%TimeScaleSlider.value)
	_refresh_row_states()

	# Toolbar buttons
	%PauseBtn.pressed.connect(_on_pause_btn_pressed)
	%GhostBtn.pressed.connect(_on_ghost_btn_pressed)
	%TraceBtn.pressed.connect(_on_trace_btn_pressed)
	%ClearBtn.pressed.connect(_on_clear_btn_pressed)
	%LinearBtn.pressed.connect(_on_linear_btn_pressed)
	%DirBtn.pressed.connect(_on_dir_btn_pressed)
	%ResetBtn.pressed.connect(_on_reset_pressed)

	if _camera_rig:
		%FrontBtn.pressed.connect(_camera_rig.snap_front)
		%SideBtn.pressed.connect(_camera_rig.snap_right)
		%TopBtn.pressed.connect(_camera_rig.snap_top)
		%IsoBtn.pressed.connect(_camera_rig.snap_iso)
		%HomeBtn.pressed.connect(_camera_rig.snap_default)
		%FitBtn.pressed.connect(_camera_rig.auto_fit)


func _build_tier_tabs() -> void:
	_tab_container = TabContainer.new()
	_tab_container.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	_tab_container.size_flags_vertical = Control.SIZE_EXPAND_FILL

	var tier_names: Array = ["Charge Photon (1-4)", "High Photon (5-8)", "Meson (9-12)", "Baryon (13-16)"]

	for tier in TIER_COUNT:
		var tab_vbox: VBoxContainer = VBoxContainer.new()
		tab_vbox.name = tier_names[tier]
		tab_vbox.add_theme_constant_override("separation", 4)

		if tier == 0:
			var help: Label = Label.new()
			help.modulate = Color(0.75, 0.75, 0.78, 1)
			help.text = "Drag any slider to activate. Prior levels auto-set to +c."
			tab_vbox.add_child(help)

		for j in LEVELS_PER_TIER:
			var level: int = tier * LEVELS_PER_TIER + j + 1

			var row: VBoxContainer = VBoxContainer.new()
			row.add_theme_constant_override("separation", 2)

			var lbl: Label = Label.new()
			lbl.text = _format_level_label(level)
			row.add_child(lbl)
			_level_labels.append(lbl)

			var slider_row: HBoxContainer = HBoxContainer.new()

			var rev_btn: Button = Button.new()
			rev_btn.text = "±"
			rev_btn.custom_minimum_size = Vector2(32, 0)
			rev_btn.tooltip_text = "Reverse spin direction"
			slider_row.add_child(rev_btn)
			_reverse_buttons.append(rev_btn)

			var slider: HSlider = HSlider.new()
			slider.min_value = -1.0
			slider.max_value = 1.0
			slider.step = 0.001
			slider.size_flags_horizontal = Control.SIZE_EXPAND_FILL
			slider.editable = true
			slider.tick_count = 5
			slider.ticks_on_borders = true
			slider_row.add_child(slider)
			_sliders.append(slider)

			var val_lbl: Label = Label.new()
			val_lbl.custom_minimum_size = Vector2(80, 0)
			val_lbl.horizontal_alignment = HORIZONTAL_ALIGNMENT_RIGHT
			val_lbl.text = "—"
			slider_row.add_child(val_lbl)
			_value_labels.append(val_lbl)

			row.add_child(slider_row)
			tab_vbox.add_child(row)

		_tab_container.add_child(tab_vbox)

	var vbox: VBoxContainer = $ControlPanel/Margin/VBox
	vbox.add_child(_tab_container)
	vbox.move_child(_tab_container, 1)


func _process(_delta: float) -> void:
	if _focus == null:
		return

	var active: int = int(_focus.call("level_count"))
	if active != _prev_level_count:
		_prev_level_count = active
		if _trace:
			_trace.clear_trace()
		if _camera_rig:
			_camera_rig.auto_fit()
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

	# Readouts
	%ParticleLabel.text = String(_focus.call("classification"))
	%SignatureLabel.text = String(_focus.call("signature"))
	%RadiusLabel.text = "%d r" % int(_focus.call("effective_radius"))
	if active > 0:
		var top_v: float = float(_focus.call("get_level_velocity", active))
		%TopOmegaLabel.text = "%+0.3f c (lvl %d)" % [top_v, active]
	else:
		%TopOmegaLabel.text = "—"
	%FpsLabel.text = "%d" % Engine.get_frames_per_second()

	# Linear state
	var linear_on: bool = bool(_focus.call("is_linear_enabled"))
	%LinearSpeedSlider.editable = linear_on
	if linear_on:
		var spd: float = float(_focus.call("get_linear_speed"))
		if not %LinearSpeedSlider.has_focus():
			%LinearSpeedSlider.set_value_no_signal(spd)
		%LinearSpeedValue.text = "%.3f c" % spd
		var dir_label: String = String(_focus.call("get_linear_direction_label"))
		%LinearVelLabel.text = "%.3f c (%s)" % [spd, dir_label]
	else:
		%LinearSpeedValue.text = "—"
		%LinearVelLabel.text = "at rest"
	var wl: float = float(_focus.call("get_wavelength"))
	if wl > 0.01:
		%WavelengthLabel.text = "%.2f r" % wl
	else:
		var char_wl: float = float(_focus.call("get_characteristic_wavelength"))
		if char_wl > 0.01:
			%WavelengthLabel.text = "%.2f r (at c)" % char_wl
		else:
			%WavelengthLabel.text = "—"

	# SI wavelength and EM band
	var wl_si: String = String(_focus.call("get_wavelength_si"))
	%WavelengthSILabel.text = wl_si
	%EMBandLabel.text = String(_focus.call("get_em_band"))
	var em_color: Color = Color(_focus.call("get_wavelength_color"))
	%EMColorSwatch.color = em_color

	# Toolbar state labels
	%PauseBtn.text = "Resume" if _paused else "Pause"
	%LinearBtn.text = "Linear: on" if linear_on else "Linear: off"
	var dir_lbl: String = String(_focus.call("get_linear_direction_label"))
	%DirBtn.text = "Dir: " + dir_lbl
	%DirBtn.disabled = not linear_on
	if _ghost:
		%GhostBtn.text = "Ghost: " + _ghost.display_mode_label()
	if _trace:
		%TraceBtn.text = "Trace: " + _trace.display_mode_label()


func _on_slider_changed(value: float, level: int) -> void:
	var snapped: float = _snap_to_targets(value, [-1.0, 0.0, 1.0], 0.02)
	if snapped != value:
		_sliders[level - 1].set_value_no_signal(snapped)
		value = snapped
	if _focus == null:
		return
	var active: int = int(_focus.call("level_count"))
	if level > active:
		for prev in range(1, level):
			if prev > active:
				_focus.call("activate_next")
			_focus.call("set_level_velocity", prev, 1.0)
			_sliders[prev - 1].set_value_no_signal(1.0)
		if int(_focus.call("level_count")) < level:
			_focus.call("activate_next")
		_refresh_row_states()
	_focus.call("set_level_velocity", level, value)


func _on_reverse_pressed(level: int) -> void:
	if _focus == null:
		return
	var active: int = int(_focus.call("level_count"))
	if level > active:
		return
	var v: float = float(_focus.call("get_level_velocity", level))
	_focus.call("set_level_velocity", level, -v)
	_sliders[level - 1].set_value_no_signal(-v)


func _on_time_scale_changed(slider_log: float) -> void:
	var snapped: float = _snap_to_targets(slider_log, [-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0], 0.02)
	if snapped != slider_log:
		%TimeScaleSlider.set_value_no_signal(snapped)
		slider_log = snapped
	_apply_time_scale(slider_log)


func _snap_to_targets(value: float, targets: Array, deadzone: float) -> float:
	for target: float in targets:
		if abs(value - target) <= deadzone:
			return target
	return value


func _on_linear_speed_changed(value: float) -> void:
	if _focus == null:
		return
	_focus.call("set_linear_speed", value)


# ---- Toolbar handlers ----

func _on_pause_btn_pressed() -> void:
	_set_paused(not _paused)


func _on_ghost_btn_pressed() -> void:
	if _ghost:
		_ghost.cycle_display_mode()


func _on_trace_btn_pressed() -> void:
	if _trace:
		_trace.cycle_display_mode()


func _on_clear_btn_pressed() -> void:
	if _trace:
		_trace.clear_trace()


func _on_linear_btn_pressed() -> void:
	if _focus:
		var current: bool = bool(_focus.call("is_linear_enabled"))
		_focus.call("set_linear_enabled", not current)


func _on_dir_btn_pressed() -> void:
	if _focus:
		var current: int = int(_focus.call("get_linear_direction_preset"))
		var count: int = int(_focus.call("get_linear_direction_count"))
		_focus.call("set_linear_direction_preset", (current + 1) % count)


func _unhandled_key_input(event: InputEvent) -> void:
	if not (event is InputEventKey) or not event.pressed or event.echo:
		return
	match (event as InputEventKey).keycode:
		KEY_SPACE:
			_set_paused(not _paused)
		KEY_L:
			if _focus != null:
				var current: bool = bool(_focus.call("is_linear_enabled"))
				_focus.call("set_linear_enabled", not current)


func _set_paused(p: bool) -> void:
	_paused = p
	if _focus != null:
		_focus.process_mode = Node.PROCESS_MODE_DISABLED if _paused else Node.PROCESS_MODE_INHERIT
	%PauseLabel.text = "PAUSED" if _paused else ""


func _on_reset_pressed() -> void:
	if _focus == null:
		return
	_focus.call("reset_stack")
	for i in PHASE_MAX_LEVEL:
		_sliders[i].set_value_no_signal(0.0)
	%TimeScaleSlider.set_value_no_signal(0.0)
	_apply_time_scale(0.0)
	%LinearSpeedSlider.set_value_no_signal(0.0)
	%LinearSpeedSlider.editable = false
	_refresh_row_states()
	_tab_container.current_tab = 0
	if _camera_rig:
		_camera_rig.snap_default()


func _apply_time_scale(slider_log: float) -> void:
	var factor: float = pow(10.0, slider_log)
	var ts: float = BASE_TIME_SCALE * factor
	if _focus != null:
		_focus.call("set_time_scale", ts)
	%TimeScaleValue.text = "%.3fx (%.2f rad/s @ c)" % [factor, ts]


func _refresh_row_states() -> void:
	if _focus == null:
		return
	var active: int = int(_focus.call("level_count"))
	for i in PHASE_MAX_LEVEL:
		var level: int = i + 1
		var is_active: bool = level <= active
		_sliders[i].modulate = Color(1, 1, 1, 1) if is_active else Color(1, 1, 1, 0.5)
		_reverse_buttons[i].disabled = not is_active

	for tier in TIER_COUNT:
		_tab_container.set_tab_disabled(tier, false)


func _format_level_label(level: int) -> String:
	if _focus == null:
		return "Lvl %d" % level
	var spin_type: String = String(_focus.call("spin_type_label", level))
	var tier: String = String(_focus.call("tier_label", level))
	var amp: float = float(_focus.call("level_amplitude", level))
	return "Lvl %d — %s (%s, r=%d)" % [level, spin_type, tier, int(amp)]
