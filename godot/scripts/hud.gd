extends CanvasLayer
## M4 HUD: toolbar + tier-tabbed slider panel for the spin stack (levels 1-15).

const PHASE_MAX_LEVEL: int = 15
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
var _field_count_spin: SpinBox
var _anticharge_slider: HSlider
var _anticharge_label: Label
var _dir_buttons: Dictionary = {}
var _exit_holes_btn: Button
var _heatmap_btn: Button
var _zone_btn: Button
var _sprinkler_btn: Button
var _capsule_btn: Button
var _decay_btn: Button
var _blank_btn: Button
var _emission_label: Label
var _preset_dropdown: OptionButton
var _presets: Array = []
var _tab_container: TabContainer
var _level_labels: Array[Label] = []
var _sliders: Array[HSlider] = []
var _value_labels: Array[Label] = []
var _spinboxes: Array[SpinBox] = []
var _reverse_buttons: Array[Button] = []
var _annul_buttons: Array[Button] = []
var _paused: bool = false
var _prev_level_count: int = 0
var _angle_graph: Control
var _angle_graph_panel: PanelContainer
var _graph_update_tick: int = 0


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
		_spinboxes[i].value_changed.connect(_on_spinbox_changed.bind(level))
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

	_build_field_panel()
	_build_angle_graph_panel()

	# Sync default button states to field_sim
	if _field_sim:
		_apply_field_direction()
		_field_sim.set_show_heatmap(true)
		_field_sim.call("set_chirality_strength", 0.15)


func _build_tier_tabs() -> void:
	_tab_container = TabContainer.new()
	_tab_container.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	_tab_container.size_flags_vertical = Control.SIZE_EXPAND_FILL

	var tier_names: Array = ["Charge Photon (1-4)", "High Photon (5-8)", "Electron/Meson (9-12)", "Uberon (13-15)"]

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
			if level > PHASE_MAX_LEVEL:
				break

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
			val_lbl.custom_minimum_size = Vector2(48, 0)
			val_lbl.horizontal_alignment = HORIZONTAL_ALIGNMENT_RIGHT
			val_lbl.text = "—"
			slider_row.add_child(val_lbl)
			_value_labels.append(val_lbl)

			var spinbox: SpinBox = SpinBox.new()
			spinbox.min_value = -1.0
			spinbox.max_value = 1.0
			spinbox.step = 0.001
			spinbox.custom_minimum_size = Vector2(80, 0)
			spinbox.allow_greater = false
			spinbox.allow_lesser = false
			spinbox.tooltip_text = "Type exact velocity (c units)"
			slider_row.add_child(spinbox)
			_spinboxes.append(spinbox)

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
			if not _spinboxes[i].get_line_edit().has_focus():
				_spinboxes[i].set_value_no_signal(v)
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
	_update_angle_graph()

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
	var snapped: float = _snap_to_targets(value, [-1.0, -0.5, 0.0, 0.5, 1.0], 0.02)
	if snapped != value:
		_sliders[level - 1].set_value_no_signal(snapped)
		value = snapped
	_spinboxes[level - 1].set_value_no_signal(value)
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


func _on_spinbox_changed(value: float, level: int) -> void:
	_sliders[level - 1].set_value_no_signal(value)
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
	_spinboxes[level - 1].set_value_no_signal(-v)


func _on_annul_pressed(level: int) -> void:
	if _focus == null:
		return
	var active: int = int(_focus.call("level_count"))
	if level > active:
		return
	_focus.call("truncate_to_level", level)
	for i in range(level - 1, PHASE_MAX_LEVEL):
		_sliders[i].set_value_no_signal(0.0)
		_spinboxes[i].set_value_no_signal(0.0)
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
		_spinboxes[i].set_value_no_signal(0.0)
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
		_spinboxes[i].editable = true
		_spinboxes[i].modulate = Color(1, 1, 1, 1) if is_active else Color(1, 1, 1, 0.5)
		_reverse_buttons[i].disabled = not is_active
		_annul_buttons[i].disabled = not is_active

	for tier in TIER_COUNT:
		_tab_container.set_tab_disabled(tier, false)


func _build_field_panel() -> void:
	var panel: PanelContainer = PanelContainer.new()
	panel.anchor_left = 0.0
	panel.anchor_right = 0.0
	panel.anchor_top = 1.0
	panel.anchor_bottom = 1.0
	panel.offset_left = 16
	panel.offset_right = 340
	panel.offset_top = -310
	panel.offset_bottom = -16
	panel.mouse_filter = Control.MOUSE_FILTER_STOP

	var margin: MarginContainer = MarginContainer.new()
	margin.add_theme_constant_override("margin_left", 10)
	margin.add_theme_constant_override("margin_top", 8)
	margin.add_theme_constant_override("margin_right", 10)
	margin.add_theme_constant_override("margin_bottom", 8)
	panel.add_child(margin)

	var vbox: VBoxContainer = VBoxContainer.new()
	vbox.add_theme_constant_override("separation", 4)
	margin.add_child(vbox)

	# Title
	var title: Label = Label.new()
	title.text = "Charge Field"
	title.add_theme_font_size_override("font_size", 15)
	vbox.add_child(title)

	var sep_top: HSeparator = HSeparator.new()
	vbox.add_child(sep_top)

	# Row 1: Start / Stop / Reset + status
	var row1: HBoxContainer = HBoxContainer.new()
	row1.add_theme_constant_override("separation", 4)

	_field_start_btn = Button.new()
	_field_start_btn.text = "Start"
	_field_start_btn.pressed.connect(_on_field_start_pressed)
	row1.add_child(_field_start_btn)

	_field_stop_btn = Button.new()
	_field_stop_btn.text = "Stop"
	_field_stop_btn.disabled = true
	_field_stop_btn.pressed.connect(_on_field_stop_pressed)
	row1.add_child(_field_stop_btn)

	_field_reset_btn = Button.new()
	_field_reset_btn.text = "Reset"
	_field_reset_btn.disabled = true
	_field_reset_btn.pressed.connect(_on_field_reset_pressed)
	row1.add_child(_field_reset_btn)

	_field_count_label = Label.new()
	_field_count_label.text = "off"
	_field_count_label.modulate = Color(0.7, 0.8, 0.9, 1)
	row1.add_child(_field_count_label)

	vbox.add_child(row1)

	# Row 2: Photon count
	var row2: HBoxContainer = HBoxContainer.new()
	row2.add_theme_constant_override("separation", 4)

	var count_lbl: Label = Label.new()
	count_lbl.text = "Photons:"
	row2.add_child(count_lbl)

	_field_count_spin = SpinBox.new()
	_field_count_spin.min_value = 50
	_field_count_spin.max_value = 100000
	_field_count_spin.step = 50
	_field_count_spin.value = 10000
	_field_count_spin.custom_minimum_size = Vector2(90, 0)
	_field_count_spin.tooltip_text = "Number of charge photons (50-100k)"
	row2.add_child(_field_count_spin)

	var ac_lbl: Label = Label.new()
	ac_lbl.text = "  AC:"
	row2.add_child(ac_lbl)

	_anticharge_slider = HSlider.new()
	_anticharge_slider.min_value = 0.0
	_anticharge_slider.max_value = 100.0
	_anticharge_slider.step = 1.0
	_anticharge_slider.value = 33.0
	_anticharge_slider.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	_anticharge_slider.custom_minimum_size = Vector2(50, 0)
	_anticharge_slider.tooltip_text = "Percentage of anti-charge photons"
	_anticharge_slider.value_changed.connect(_on_anticharge_changed)
	row2.add_child(_anticharge_slider)

	_anticharge_label = Label.new()
	_anticharge_label.custom_minimum_size = Vector2(34, 0)
	_anticharge_label.text = "33%"
	row2.add_child(_anticharge_label)

	vbox.add_child(row2)

	# Direction row
	var dir_row: HBoxContainer = HBoxContainer.new()
	dir_row.add_theme_constant_override("separation", 2)

	var dir_lbl: Label = Label.new()
	dir_lbl.text = "Dir:"
	dir_row.add_child(dir_lbl)

	var axes: Array = ["+X", "-X", "+Y", "-Y", "+Z", "-Z"]
	var default_dirs: Array = ["+Z", "-Z"]
	for axis_name: String in axes:
		var btn: Button = Button.new()
		btn.text = axis_name
		btn.toggle_mode = true
		btn.button_pressed = axis_name in default_dirs
		btn.custom_minimum_size = Vector2(32, 0)
		btn.pressed.connect(_on_direction_toggled)
		dir_row.add_child(btn)
		_dir_buttons[axis_name] = btn

	vbox.add_child(dir_row)

	var sep_mid: HSeparator = HSeparator.new()
	vbox.add_child(sep_mid)

	# Toggle grid: 2 columns of CheckButtons
	var toggle_grid: GridContainer = GridContainer.new()
	toggle_grid.columns = 2
	toggle_grid.add_theme_constant_override("h_separation", 16)
	toggle_grid.add_theme_constant_override("v_separation", 2)

	_heatmap_btn = CheckButton.new()
	_heatmap_btn.text = "Heatmap"
	_heatmap_btn.button_pressed = true
	_heatmap_btn.tooltip_text = "Accumulate emission pattern as heatmap on sphere"
	_heatmap_btn.pressed.connect(_on_heatmap_toggled)
	toggle_grid.add_child(_heatmap_btn)

	_zone_btn = CheckButton.new()
	_zone_btn.text = "Zone"
	_zone_btn.button_pressed = true
	_zone_btn.tooltip_text = "Zone collision: polar pass-through + equatorial deflection"
	_zone_btn.pressed.connect(_on_zone_toggled)
	toggle_grid.add_child(_zone_btn)

	_sprinkler_btn = CheckButton.new()
	_sprinkler_btn.text = "Emission"
	_sprinkler_btn.button_pressed = true
	_sprinkler_btn.tooltip_text = "Pole emission: charge enters at spin poles, exits equator"
	_sprinkler_btn.pressed.connect(_on_sprinkler_toggled)
	toggle_grid.add_child(_sprinkler_btn)

	_decay_btn = CheckButton.new()
	_decay_btn.text = "Decay"
	_decay_btn.button_pressed = false
	_decay_btn.tooltip_text = "Fade old heatmap hits over time (off = accumulate forever)"
	_decay_btn.pressed.connect(_on_decay_toggled)
	toggle_grid.add_child(_decay_btn)

	_capsule_btn = CheckButton.new()
	_capsule_btn.text = "Capsules"
	_capsule_btn.button_pressed = false
	_capsule_btn.tooltip_text = "Show trace capsule collision geometry"
	_capsule_btn.pressed.connect(_on_capsule_toggled)
	toggle_grid.add_child(_capsule_btn)

	_blank_btn = CheckButton.new()
	_blank_btn.text = "Blank"
	_blank_btn.button_pressed = false
	_blank_btn.tooltip_text = "Hide charge cloud (exit dots remain)"
	_blank_btn.pressed.connect(_on_blank_toggled)
	toggle_grid.add_child(_blank_btn)

	_exit_holes_btn = CheckButton.new()
	_exit_holes_btn.text = "Exit dots"
	_exit_holes_btn.button_pressed = false
	_exit_holes_btn.tooltip_text = "Show exit points on ghost sphere"
	_exit_holes_btn.pressed.connect(_on_exit_holes_toggled)
	toggle_grid.add_child(_exit_holes_btn)

	vbox.add_child(toggle_grid)

	_emission_label = Label.new()
	_emission_label.text = ""
	_emission_label.modulate = Color(0.8, 0.9, 0.7, 1)
	vbox.add_child(_emission_label)

	add_child(panel)


func _build_angle_graph_panel() -> void:
	_angle_graph_panel = PanelContainer.new()
	_angle_graph_panel.anchor_left = 1.0
	_angle_graph_panel.anchor_right = 1.0
	_angle_graph_panel.anchor_top = 1.0
	_angle_graph_panel.anchor_bottom = 1.0
	_angle_graph_panel.offset_left = -340
	_angle_graph_panel.offset_right = -16
	_angle_graph_panel.offset_top = -180
	_angle_graph_panel.offset_bottom = -16
	_angle_graph_panel.mouse_filter = Control.MOUSE_FILTER_STOP
	_angle_graph_panel.visible = false

	var graph_script: GDScript = load("res://scripts/exit_angle_graph.gd")
	_angle_graph = Control.new()
	_angle_graph.set_script(graph_script)
	_angle_graph.set_anchors_preset(Control.PRESET_FULL_RECT)
	_angle_graph_panel.add_child(_angle_graph)

	add_child(_angle_graph_panel)


func _on_field_start_pressed() -> void:
	if _field_sim:
		if _focus:
			var outermost: int = int(_focus.call("get_outermost_level"))
			if outermost < 9:
				_emission_label.text = "Field requires baryon tier (level 9+)"
				return
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
		_field_count_label.text = "off"
	_reset_field_controls()


func _reset_field_controls() -> void:
	_field_start_btn.disabled = false
	_field_stop_btn.disabled = true
	_field_reset_btn.disabled = true
	_field_count_label.text = "off"
	var default_dirs: Array = ["+Z", "-Z"]
	for key: String in _dir_buttons:
		_dir_buttons[key].button_pressed = key in default_dirs
	_exit_holes_btn.button_pressed = false
	_heatmap_btn.button_pressed = true
	_blank_btn.button_pressed = false
	_zone_btn.button_pressed = true
	_sprinkler_btn.button_pressed = true
	_capsule_btn.button_pressed = false
	_decay_btn.button_pressed = false
	_anticharge_slider.set_value_no_signal(33.0)
	_anticharge_label.text = "33%"
	if _field_sim:
		_field_sim.set_direction_strength(0.0)
		_field_sim.set_direction_scatter(0.0)
		_field_sim.set_photon_ratio(0.67)
		_field_sim.set_show_exit_holes(false)
		_field_sim.set_show_heatmap(true)
		_field_sim.set_heatmap_decay(false)
		_field_sim.set_cloud_blanked(false)
		_field_sim.set_zone_collision(true)
		_field_sim.set_sprinkler_emission(true)
		_field_sim.set_show_capsules(false)
		_apply_field_direction()


func _on_anticharge_changed(value: float) -> void:
	var pct: int = int(value)
	_anticharge_label.text = "%d%%" % pct
	if _field_sim:
		_field_sim.set_photon_ratio(1.0 - value / 100.0)


func _on_direction_toggled() -> void:
	_apply_field_direction()


func _on_exit_holes_toggled() -> void:
	if _field_sim:
		_field_sim.set_show_exit_holes(_exit_holes_btn.button_pressed)


func _on_blank_toggled() -> void:
	if _field_sim:
		_field_sim.set_cloud_blanked(_blank_btn.button_pressed)


func _on_heatmap_toggled() -> void:
	if _field_sim:
		_field_sim.set_show_heatmap(_heatmap_btn.button_pressed)




func _on_zone_toggled() -> void:
	if _field_sim:
		_field_sim.set_zone_collision(_zone_btn.button_pressed)


func _on_sprinkler_toggled() -> void:
	if _field_sim:
		_field_sim.set_sprinkler_emission(_sprinkler_btn.button_pressed)


func _on_capsule_toggled() -> void:
	if _field_sim:
		_field_sim.set_show_capsules(_capsule_btn.button_pressed)


func _on_decay_toggled() -> void:
	if _field_sim:
		_field_sim.set_heatmap_decay(_decay_btn.button_pressed)


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
	else:
		_field_sim.set_direction_strength(1.0)
	_field_sim.set_direction_scatter(0.0)


func _show_heatmap_readout() -> bool:
	return _heatmap_btn.button_pressed or _exit_holes_btn.button_pressed


func _update_field_readout() -> void:
	if not _field_sim:
		return

	if _field_sim.is_running():
		var count: int = _field_sim.get_photon_count()
		_field_count_label.text = "%dK" % (count / 1000)
		if _show_heatmap_readout():
			var angle: float = _field_sim.get_mean_exit_latitude()
			var hits: int = _field_sim.get_heatmap_total_hits()
			_emission_label.text = "Mean lat: %.1f° (n=%d)" % [angle, hits]
		else:
			_emission_label.text = ""
	else:
		_field_count_label.text = "off"
		_emission_label.text = ""


func _update_angle_graph() -> void:
	if _angle_graph == null or _angle_graph_panel == null:
		return
	if not _field_sim:
		_angle_graph_panel.visible = false
		return

	var show: bool = _field_sim.is_running() and _show_heatmap_readout()
	_angle_graph_panel.visible = show
	if not show:
		return

	# Update every 6th frame (~10 Hz) to save compute
	_graph_update_tick += 1
	if _graph_update_tick % 6 != 0:
		return

	var histogram: PackedFloat32Array = _field_sim.get_exit_angle_histogram()
	var sample_n: int = _field_sim.get_exit_angle_sample_count()
	var mean_deg: float = _field_sim.get_mean_exit_latitude()
	_angle_graph.update_data(histogram, sample_n, mean_deg)


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
		_spinboxes[i].set_value_no_signal(0.0)

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
		_spinboxes[i].set_value_no_signal(v)

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
