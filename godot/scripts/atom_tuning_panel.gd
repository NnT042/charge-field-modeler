extends CanvasLayer
#
# DEBUG tuning UI — labeled sliders for the atom-mode sandbox experiment
# knobs. Hidden by default. Toggled by the HUD's Debug button.
#
# The nine LOCKED force couplings (G_q, C_q, V_q, D_q, T_q, I_q, Corot, S_q,
# P_amb) are intentionally NOT exposed here any more — they are pinned by
# the M5 test suite and derived per Couplings::default() in
# rust/src/atom_core.rs ("not sliders any more"). Fiddling them from the UI
# only breaks calibrated physics.
#
# What remains are the live experiment knobs actively used this session:
# nuclear binding v2 (channeling / nuclear ambient / intra boost) and skin
# v4 mote parameters (fraction / scale / lifetime). Contains NO physics —
# it is a pure getter/setter bridge to AtomSim.
#

var _sim: Node = null
var _mode: Node = null
var _rows: Array = []  # each: {slider, value_label, getter, setter, updating}

const PARAMS := [
	{"section": "Nuclear binding (sandbox)"},
	{"label": "Channeling",    "get": "get_channeling_coupling", "set": "set_channeling_coupling", "min": 0.0, "max": 1.0,  "step": 0.005},
	{"label": "Nuclear ambient", "get": "get_nuclear_ambient",   "set": "set_nuclear_ambient",      "min": 0.0, "max": 30.0, "step": 0.1},
	{"label": "Intra boost",   "get": "get_intra_boost",         "set": "set_intra_boost",          "min": 0.0, "max": 6.0,  "step": 0.05},
	{"section": "Skin motes"},
	{"label": "Mote fraction", "get": "get_mote_fraction",       "set": "set_mote_fraction",         "min": 0.0, "max": 1.0,  "step": 0.01},
	{"label": "Mote scale",    "get": "get_mote_scale",          "set": "set_mote_scale",            "min": 0.5, "max": 8.0,  "step": 0.1},
	{"label": "Mote lifetime", "get": "get_mote_lifetime",       "set": "set_mote_lifetime",         "min": 0.2, "max": 8.0,  "step": 0.1},
]

func setup(sim: Node, mode: Node = null) -> void:
	_sim = sim
	_mode = mode
	layer = 10

	var panel := PanelContainer.new()
	panel.custom_minimum_size = Vector2(320, 0)
	add_child(panel)

	var vbox := VBoxContainer.new()
	vbox.add_theme_constant_override("separation", 4)
	panel.add_child(vbox)

	var title := Label.new()
	title.text = "Force Tuning"
	title.add_theme_font_size_override("font_size", 13)
	vbox.add_child(title)

	for p in PARAMS:
		if p.has("section"):
			vbox.add_child(HSeparator.new())
			var section_label := Label.new()
			section_label.text = p["section"]
			section_label.add_theme_font_size_override("font_size", 12)
			vbox.add_child(section_label)
		else:
			_add_row(vbox, p)

	# Debug visualization toggles (superseded visuals, kept for inspection)
	if _mode != null:
		vbox.add_child(HSeparator.new())
		_add_debug_check(vbox, "Pole axis lines", "show_pole_lines")
		_add_debug_check(vbox, "Profile rings (wireframe)", "show_profile_rings")

	# Position: mid-right edge, biased slightly below center (0.55) so it
	# clears the top-right ReadoutPanel (~y=52..300). All rows/sections are
	# built above this point so reset_size() picks up the real combined
	# minimum height before we anchor+offset the rect — computing offsets
	# off a not-yet-laid-out size (still (0,0) at this point otherwise)
	# would collapse the panel to zero height.
	panel.reset_size()
	panel.set_anchors_and_offsets_preset(Control.PRESET_CENTER_RIGHT)
	panel.anchor_top = 0.55
	panel.anchor_bottom = 0.55
	panel.offset_right = -12.0
	panel.offset_left = panel.offset_right - 320.0

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
