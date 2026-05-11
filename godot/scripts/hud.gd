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
@export var field_sim_path: NodePath
@export var field_renderer_path: NodePath

var _focus: Node3D
var _camera_rig: Node3D
var _ghost: Node3D
var _trace: Node3D
var _field_sim: Node3D
var _field_renderer: Node3D
var _field_start_btn: Button
var _field_stop_btn: Button
var _field_reset_btn: Button
var _field_count_label: Label
var _collision_label: Label
var _field_count_spin: SpinBox
var _anticharge_slider: HSlider
var _anticharge_label: Label
var _dir_buttons: Dictionary = {}
var _parallel_btn: Button
var _exit_holes_btn: Button
var _blank_btn: Button
var _stream_btn: Button
var _preset_dropdown: OptionButton
var _presets: Array = []
var _tab_container: TabContainer
var _level_labels: Array[Label] = []
var _sliders: Array[HSlider] = []
var _value_labels: Array[Label] = []
var _reverse_buttons: Array[Button] = []
var _annul_buttons: Array[Button] = []
var _paused: bool = false
var _prev_level_count: int = 0


func _ready() -> void:
	_focus = get_node_or_null(focus_particle_path) as Node3D
	_camera_rig = get_node_or_null(camera_rig_path) as Node3D
	_ghost = get_node_or_null(ghost_sphere_path) as Node3D
	_trace = get_node_or_null(path_trace_path) as Node3D
	_field_sim = get_node_or_null(field_sim_path) as Node3D
	_field_renderer = get_node_or_null(field_renderer_path) as Node3D

	if _focus == null:
		push_error("HUD: focus_particle_path does not resolve to a Node3D")
		return

	_build_tier_tabs()
	_build_preset_row()

	for i in PHASE_MAX_LEVEL:
		var level: int = i + 1
		_sliders[i].value_changed.connect(_on_slider_changed.bind(level))
		_reverse_buttons[i].pressed.connect(_on_reverse_pressed.bind(level))
		_annul_buttons[i].pressed.connect(_on_annul_pressed.bind(level))

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

	var field_row: VBoxContainer = _build_field_controls()
	var vbox_toolbar: VBoxContainer = $ControlPanel/Margin/VBox
	vbox_toolbar.add_child(field_row)


func _build_tier_tabs() -> void:
	_tab_container = TabContainer.new()
	_tab_container.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	_tab_container.size_flags_vertical = Control.SIZE_EXPAND_FILL

	var tier_names: Array = ["Charge Photon (1-4)", "High Photon (5-8)", "Electron/Meson (9-12)", "Baryon (13-16)"]

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

			var annul_btn: Button = Button.new()
			annul_btn.text = "0"
			annul_btn.custom_minimum_size = Vector2(24, 0)
			annul_btn.tooltip_text = "Annul this spin and all above"
			slider_row.add_child(annul_btn)
			_annul_buttons.append(annul_btn)

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

	_update_field_readout()

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
		if _field_sim:
			_field_sim.clear_exit_holes()
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


func _on_annul_pressed(level: int) -> void:
	if _focus == null:
		return
	var active: int = int(_focus.call("level_count"))
	if level > active:
		return
	_focus.call("truncate_to_level", level)
	for i in range(level - 1, PHASE_MAX_LEVEL):
		_sliders[i].set_value_no_signal(0.0)
	var target_level: int = max(1, level - 1)
	var target_tier: int = (target_level - 1) / LEVELS_PER_TIER
	_tab_container.current_tab = target_tier
	_refresh_row_states()
	if _camera_rig:
		_camera_rig.auto_fit()


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
		KEY_F:
			if _field_renderer:
				_field_renderer.toggle_visible()


func _set_paused(p: bool) -> void:
	_paused = p
	if _focus != null:
		_focus.process_mode = Node.PROCESS_MODE_DISABLED if _paused else Node.PROCESS_MODE_INHERIT
	%PauseLabel.text = "PAUSED" if _paused else ""


func _on_reset_pressed() -> void:
	if _paused:
		_set_paused(false)
	if _focus == null:
		return
	_focus.call("reset_stack")
	for i in PHASE_MAX_LEVEL:
		_sliders[i].set_value_no_signal(0.0)
	%TimeScaleSlider.set_value_no_signal(0.0)
	_apply_time_scale(0.0)
	%LinearSpeedSlider.set_value_no_signal(0.0)
	%LinearSpeedSlider.editable = false
	if _field_sim:
		_field_sim.reset_field()
	_reset_field_controls()
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
		_annul_buttons[i].disabled = not is_active

	for tier in TIER_COUNT:
		_tab_container.set_tab_disabled(tier, false)


func _build_field_controls() -> VBoxContainer:
	var wrapper: VBoxContainer = VBoxContainer.new()
	wrapper.add_theme_constant_override("separation", 4)

	# Row 1: Start / Stop / Reset + status
	var row1: HBoxContainer = HBoxContainer.new()
	row1.add_theme_constant_override("separation", 6)

	_field_start_btn = Button.new()
	_field_start_btn.text = "Start Field"
	_field_start_btn.pressed.connect(_on_field_start_pressed)
	row1.add_child(_field_start_btn)

	_field_stop_btn = Button.new()
	_field_stop_btn.text = "Stop Field"
	_field_stop_btn.disabled = true
	_field_stop_btn.pressed.connect(_on_field_stop_pressed)
	row1.add_child(_field_stop_btn)

	_field_reset_btn = Button.new()
	_field_reset_btn.text = "Reset Field"
	_field_reset_btn.disabled = true
	_field_reset_btn.pressed.connect(_on_field_reset_pressed)
	row1.add_child(_field_reset_btn)

	_field_count_label = Label.new()
	_field_count_label.text = "Field: off"
	row1.add_child(_field_count_label)

	_collision_label = Label.new()
	_collision_label.text = ""
	_collision_label.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	_collision_label.horizontal_alignment = HORIZONTAL_ALIGNMENT_RIGHT
	row1.add_child(_collision_label)

	wrapper.add_child(row1)

	# Row 2: Count + Anticharge
	var row2: HBoxContainer = HBoxContainer.new()
	row2.add_theme_constant_override("separation", 6)

	var count_lbl: Label = Label.new()
	count_lbl.text = "Count:"
	row2.add_child(count_lbl)

	_field_count_spin = SpinBox.new()
	_field_count_spin.min_value = 50
	_field_count_spin.max_value = 100000
	_field_count_spin.step = 50
	_field_count_spin.value = 5000
	_field_count_spin.custom_minimum_size = Vector2(100, 0)
	_field_count_spin.tooltip_text = "Number of charge photons (50-100k)"
	row2.add_child(_field_count_spin)

	var sep: VSeparator = VSeparator.new()
	row2.add_child(sep)

	var ac_lbl: Label = Label.new()
	ac_lbl.text = "Anticharge:"
	row2.add_child(ac_lbl)

	_anticharge_slider = HSlider.new()
	_anticharge_slider.min_value = 0.0
	_anticharge_slider.max_value = 100.0
	_anticharge_slider.step = 1.0
	_anticharge_slider.value = 33.0
	_anticharge_slider.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	_anticharge_slider.custom_minimum_size = Vector2(80, 0)
	_anticharge_slider.tooltip_text = "Percentage of red (anti-charge) photons"
	_anticharge_slider.value_changed.connect(_on_anticharge_changed)
	row2.add_child(_anticharge_slider)

	_anticharge_label = Label.new()
	_anticharge_label.custom_minimum_size = Vector2(40, 0)
	_anticharge_label.text = "33%"
	row2.add_child(_anticharge_label)

	wrapper.add_child(row2)

	# Row 3: Direction toggles + Parallel switch
	var row3: HBoxContainer = HBoxContainer.new()
	row3.add_theme_constant_override("separation", 4)

	var dir_lbl: Label = Label.new()
	dir_lbl.text = "Direction:"
	row3.add_child(dir_lbl)

	var axes: Array = ["+X", "-X", "+Y", "-Y", "+Z", "-Z"]
	for axis_name: String in axes:
		var btn: Button = Button.new()
		btn.text = axis_name
		btn.toggle_mode = true
		btn.custom_minimum_size = Vector2(36, 0)
		btn.pressed.connect(_on_direction_toggled)
		row3.add_child(btn)
		_dir_buttons[axis_name] = btn

	var sep2: VSeparator = VSeparator.new()
	row3.add_child(sep2)

	_parallel_btn = Button.new()
	_parallel_btn.text = "Scatter"
	_parallel_btn.toggle_mode = true
	_parallel_btn.button_pressed = true
	_parallel_btn.tooltip_text = "ON = mostly-parallel (some scatter). OFF = all-parallel (exact direction)"
	_parallel_btn.pressed.connect(_on_parallel_toggled)
	row3.add_child(_parallel_btn)

	var sep3: VSeparator = VSeparator.new()
	row3.add_child(sep3)

	_exit_holes_btn = Button.new()
	_exit_holes_btn.text = "Exit dots"
	_exit_holes_btn.toggle_mode = true
	_exit_holes_btn.button_pressed = false
	_exit_holes_btn.tooltip_text = "Show exit points on ghost sphere"
	_exit_holes_btn.pressed.connect(_on_exit_holes_toggled)
	row3.add_child(_exit_holes_btn)

	_blank_btn = Button.new()
	_blank_btn.text = "Blank"
	_blank_btn.toggle_mode = true
	_blank_btn.button_pressed = false
	_blank_btn.tooltip_text = "Hide charge cloud (exit dots remain)"
	_blank_btn.pressed.connect(_on_blank_toggled)
	row3.add_child(_blank_btn)

	_stream_btn = Button.new()
	_stream_btn.text = "Stream"
	_stream_btn.toggle_mode = true
	_stream_btn.button_pressed = false
	_stream_btn.tooltip_text = "Continuous flow instead of synchronized waves"
	_stream_btn.pressed.connect(_on_stream_toggled)
	row3.add_child(_stream_btn)

	wrapper.add_child(row3)
	return wrapper


func _on_field_start_pressed() -> void:
	if _field_sim:
		var count: int = int(_field_count_spin.value)
		_field_sim.start_field(count)
		_field_start_btn.disabled = true
		_field_stop_btn.disabled = false
		_field_reset_btn.disabled = false


func _on_field_stop_pressed() -> void:
	if _field_sim:
		_field_sim.stop_field()
		_field_start_btn.disabled = false
		_field_stop_btn.disabled = true


func _on_field_reset_pressed() -> void:
	if _field_sim:
		_field_sim.reset_field()
		_field_start_btn.disabled = false
		_field_stop_btn.disabled = true
		_field_reset_btn.disabled = true
		_field_count_label.text = "Field: off"
	_reset_field_controls()


func _reset_field_controls() -> void:
	for key: String in _dir_buttons:
		_dir_buttons[key].button_pressed = false
	_parallel_btn.button_pressed = true
	_parallel_btn.text = "Scatter"
	_exit_holes_btn.button_pressed = false
	_blank_btn.button_pressed = false
	_stream_btn.button_pressed = false
	_anticharge_slider.set_value_no_signal(33.0)
	_anticharge_label.text = "33%"
	if _field_sim:
		_field_sim.set_direction_mask(0)
		_field_sim.set_direction_strength(0.0)
		_field_sim.set_direction_scatter(0.3)
		_field_sim.set_photon_ratio(0.67)
		_field_sim.set_show_exit_holes(false)
		_field_sim.set_cloud_blanked(false)
		_field_sim.set_stream_mode(false)


func _on_anticharge_changed(value: float) -> void:
	var pct: int = int(value)
	_anticharge_label.text = "%d%%" % pct
	if _field_sim:
		_field_sim.set_photon_ratio(1.0 - value / 100.0)


func _on_direction_toggled() -> void:
	_apply_field_direction()


func _on_parallel_toggled() -> void:
	_parallel_btn.text = "Scatter" if _parallel_btn.button_pressed else "Parallel"
	_apply_field_direction()


func _on_exit_holes_toggled() -> void:
	if _field_sim:
		_field_sim.set_show_exit_holes(_exit_holes_btn.button_pressed)


func _on_blank_toggled() -> void:
	if _field_sim:
		_field_sim.set_cloud_blanked(_blank_btn.button_pressed)


func _on_stream_toggled() -> void:
	if _field_sim:
		_field_sim.set_stream_mode(_stream_btn.button_pressed)


func _apply_field_direction() -> void:
	if not _field_sim:
		return

	# Bitmask: +X=0, -X=1, +Y=2, -Y=3, +Z=4, -Z=5
	var mask: int = 0
	var axis_bits: Dictionary = {"+X": 0, "-X": 1, "+Y": 2, "-Y": 3, "+Z": 4, "-Z": 5}
	for key: String in axis_bits:
		if _dir_buttons[key].button_pressed:
			mask |= 1 << axis_bits[key]

	_field_sim.set_direction_mask(mask)

	if mask == 0:
		_field_sim.set_direction_strength(0.0)
		_field_sim.set_direction_scatter(0.3)
	else:
		_field_sim.set_direction_strength(1.0)
		if _parallel_btn.button_pressed:
			_field_sim.set_direction_scatter(0.3)
		else:
			_field_sim.set_direction_scatter(0.0)


func _update_field_readout() -> void:
	if not _field_sim:
		return

	if _field_sim.is_running():
		var count: int = _field_sim.get_photon_count()
		_field_count_label.text = "Field: %dK" % (count / 1000)
		var collisions: int = _field_sim.get_collision_count()
		var impulse: Vector3 = _field_sim.get_net_impulse()
		_collision_label.text = "Cloud hits: %d  |F|: %.3f" % [collisions, impulse.length()]
	else:
		_field_count_label.text = "Field: off"
		_collision_label.text = ""


func _format_level_label(level: int) -> String:
	if _focus == null:
		return "Lvl %d" % level
	var spin_type: String = String(_focus.call("spin_type_label", level))
	var tier: String = String(_focus.call("tier_label", level))
	var amp: float = float(_focus.call("level_amplitude", level))
	return "Lvl %d — %s (%s, r=%d)" % [level, spin_type, tier, int(amp)]


# ---- Particle presets ----

func _build_preset_row() -> void:
	_presets = _load_presets()
	if _presets.is_empty():
		return

	var row: HBoxContainer = HBoxContainer.new()
	row.add_theme_constant_override("separation", 6)

	var lbl: Label = Label.new()
	lbl.text = "Preset:"
	row.add_child(lbl)

	_preset_dropdown = OptionButton.new()
	_preset_dropdown.add_item("Custom")
	for p: Dictionary in _presets:
		_preset_dropdown.add_item(p.get("name", "?"))
	_preset_dropdown.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	_preset_dropdown.item_selected.connect(_on_preset_selected)
	row.add_child(_preset_dropdown)

	var vbox: VBoxContainer = $ControlPanel/Margin/VBox
	vbox.add_child(row)
	vbox.move_child(row, 1)


func _load_presets() -> Array:
	var f: FileAccess = FileAccess.open("res://config/particle_presets.json", FileAccess.READ)
	if f == null:
		push_warning("No particle_presets.json found")
		return []
	var text: String = f.get_as_text()
	f.close()
	var json: JSON = JSON.new()
	var err: int = json.parse(text)
	if err != OK:
		push_error("particle_presets.json parse error: " + json.get_error_message())
		return []
	if json.data is Array:
		return json.data
	return []


func _on_preset_selected(index: int) -> void:
	if index == 0 or _focus == null:
		return
	var preset: Dictionary = _presets[index - 1]
	var levels: Array = preset.get("levels", [])
	if levels.is_empty():
		return

	_focus.call("reset_stack")
	for i in PHASE_MAX_LEVEL:
		_sliders[i].set_value_no_signal(0.0)

	var count: int = mini(levels.size(), PHASE_MAX_LEVEL)
	for i in count:
		var v: float = float(levels[i])
		var level: int = i + 1
		var active: int = int(_focus.call("level_count"))
		if level > active:
			var top: int = active
			if top >= 1:
				var top_v: float = float(_focus.call("get_level_velocity", top))
				if absf(top_v) < 1.0 - 1e-6:
					_focus.call("set_level_velocity", top, signf(top_v) if absf(top_v) > 0.01 else 1.0)
			_focus.call("activate_next")
		_focus.call("set_level_velocity", level, v)
		_sliders[i].set_value_no_signal(v)

	if preset.has("time_scale_log"):
		var ts_log: float = float(preset["time_scale_log"])
		%TimeScaleSlider.set_value_no_signal(ts_log)
		_apply_time_scale(ts_log)

	var target_tier: int = (count - 1) / LEVELS_PER_TIER
	_tab_container.current_tab = target_tier
	_refresh_row_states()
	if _field_sim:
		_field_sim.clear_exit_holes()
	if _trace:
		_trace.clear_trace()
	if _camera_rig:
		_camera_rig.auto_fit()
