extends CanvasLayer
#
# DEBUG tuning UI — labeled sliders for the atom-mode force couplings.
# Hidden by default (constants are locked and derived — see
# Couplings::default() in rust/src/atom_core.rs). Toggled by the HUD's
# Debug button for physics experiments. Contains NO physics — it is a pure
# getter/setter bridge to AtomSim.
#

var _sim: Node = null
var _mode: Node = null
var _rows: Array = []  # each: {slider, value_label, getter, setter, updating}

const PARAMS := [
	{"label": "G_q  Gravity",  "get": "get_gravity_coupling", "set": "set_gravity_coupling", "min": 0.0, "max": 3.0,  "step": 0.005},
	{"label": "C_q  Charge",   "get": "get_charge_coupling",  "set": "set_charge_coupling",  "min": 0.0, "max": 2000.0, "step": 5.0},
	{"label": "V_q  Vortex",   "get": "get_vortex_coupling",  "set": "set_vortex_coupling",  "min": 0.0, "max": 2.0,  "step": 0.005},
	{"label": "D_q  Doppler",  "get": "get_drag_coupling",    "set": "set_drag_coupling",    "min": 0.0, "max": 2.0,  "step": 0.005},
	{"label": "T_q  Torque",   "get": "get_torque_coupling",  "set": "set_torque_coupling",  "min": 0.0, "max": 2.0,  "step": 0.005},
	{"label": "I_q  Intake",   "get": "get_intake_coupling",  "set": "set_intake_coupling",  "min": 0.0, "max": 3.0,  "step": 0.005},
	{"label": "Corot Drag",    "get": "get_corot_coupling",   "set": "set_corot_coupling",   "min": 0.0, "max": 3.0,  "step": 0.005},
	{"label": "S_q  Stream",   "get": "get_stream_coupling",  "set": "set_stream_coupling",  "min": 0.0, "max": 100.0, "step": 0.5},
	{"label": "P_amb Ambient", "get": "get_ambient_pressure", "set": "set_ambient_pressure", "min": 0.0, "max": 1.0,  "step": 0.002},
]

func setup(sim: Node, mode: Node = null) -> void:
	_sim = sim
	_mode = mode
	layer = 10

	var panel := PanelContainer.new()
	panel.position = Vector2(12, 330)
	panel.custom_minimum_size = Vector2(320, 0)
	add_child(panel)

	var vbox := VBoxContainer.new()
	vbox.add_theme_constant_override("separation", 4)
	panel.add_child(vbox)

	var title := Label.new()
	title.text = "Force Tuning"
	title.add_theme_font_size_override("font_size", 13)
	vbox.add_child(title)

	var calib_btn := Button.new()
	calib_btn.text = "Calibrate G_q  (polar orbit)"
	calib_btn.pressed.connect(_on_calibrate)
	vbox.add_child(calib_btn)

	for p in PARAMS:
		_add_row(vbox, p)

	# Debug visualization toggles (superseded visuals, kept for inspection)
	if _mode != null:
		vbox.add_child(HSeparator.new())
		_add_debug_check(vbox, "Pole axis lines", "show_pole_lines")
		_add_debug_check(vbox, "Profile rings (wireframe)", "show_profile_rings")

func _add_debug_check(parent: VBoxContainer, label: String, prop: String) -> void:
	var check := CheckBox.new()
	check.text = label
	check.add_theme_font_size_override("font_size", 12)
	check.button_pressed = _mode.get(prop)
	check.toggled.connect(func(on: bool): _mode.set(prop, on))
	parent.add_child(check)

func _add_row(parent: VBoxContainer, p: Dictionary) -> void:
	var row := HBoxContainer.new()
	row.add_theme_constant_override("separation", 6)
	parent.add_child(row)

	var name_label := Label.new()
	name_label.text = p["label"]
	name_label.custom_minimum_size = Vector2(110, 0)
	name_label.add_theme_font_size_override("font_size", 12)
	row.add_child(name_label)

	var slider := HSlider.new()
	slider.min_value = p["min"]
	slider.max_value = p["max"]
	slider.step = p["step"]
	slider.custom_minimum_size = Vector2(140, 0)
	slider.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	row.add_child(slider)

	var value_label := Label.new()
	value_label.custom_minimum_size = Vector2(48, 0)
	value_label.add_theme_font_size_override("font_size", 12)
	row.add_child(value_label)

	var entry := {"slider": slider, "value_label": value_label,
		"getter": p["get"], "setter": p["set"], "updating": false}
	slider.value_changed.connect(_on_slider_changed.bind(entry))
	_rows.append(entry)

func _on_slider_changed(value: float, entry: Dictionary) -> void:
	if entry["updating"]:
		return  # programmatic refresh, not a user drag — don't echo back
	_sim.call(entry["setter"], value)
	entry["value_label"].text = "%.3f" % value

func _on_calibrate() -> void:
	if _sim == null:
		return
	var g: float = _sim.call("auto_calibrate_polar")
	print("[tuning] auto_calibrate_polar -> G_q = %.4f" % g)

func _process(_delta: float) -> void:
	# Mirror live values (keyboard +/-, calibrate) back into the sliders,
	# but never fight a slider the user is currently dragging.
	if _sim == null:
		return
	for entry in _rows:
		var slider: HSlider = entry["slider"]
		if slider.has_focus():
			continue
		var v: float = _sim.call(entry["getter"])
		entry["updating"] = true
		slider.value = v
		entry["updating"] = false
		entry["value_label"].text = "%.3f" % v
