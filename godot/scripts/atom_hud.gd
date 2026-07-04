extends CanvasLayer
#
# Atom-mode HUD: toolbar, scenario panel, readout panel.
# Pure UI — all physics lives in Rust (AtomSim); all sim wiring lives in
# atom_mode.gd. This script only calls public methods on the mode node.
#

@export var atom_mode_path: NodePath = NodePath("..")
@export var camera_rig_path: NodePath = NodePath("../OrbitCameraRig")

var _mode: Node = null
var _camera_rig: Node3D = null

func _ready() -> void:
	_mode = get_node_or_null(atom_mode_path)
	_camera_rig = get_node_or_null(camera_rig_path)

	%PauseBtn.pressed.connect(func(): _mode.toggle_pause())
	%ResetBtn.pressed.connect(func(): _mode.reset_scenario())
	%CloudsBtn.toggled.connect(func(on: bool): _mode.set_clouds(on))
	%DebugBtn.toggled.connect(func(on: bool): _mode.set_tuning_visible(on))
	%ModeBtn.pressed.connect(func():
		get_tree().change_scene_to_file("res://scenes/main.tscn"))

	if _camera_rig:
		%FrontBtn.pressed.connect(_camera_rig.snap_front)
		%SideBtn.pressed.connect(_camera_rig.snap_right)
		%TopBtn.pressed.connect(_camera_rig.snap_top)
		%IsoBtn.pressed.connect(_camera_rig.snap_iso)

	var scenarios := {
		"ProtonsBtn": "protons",
		"HydrogenBtn": "hydrogen",
		"H2BondBtn": "h2_bond",
		"H2RepelBtn": "h2_repel",
		"AlphaBtn": "alpha",
		"HeliumBtn": "helium",
		"CarbonBtn": "carbon",
		"NitrogenBtn": "nitrogen",
		"OxygenBtn": "oxygen",
		"NeonBtn": "neon",
		"ArgonBtn": "argon",
	}
	for btn_name in scenarios:
		var node := get_node("%" + btn_name) as Button
		node.pressed.connect(_mode.spawn_scenario.bind(scenarios[btn_name]))

	%SpawnPBtn.pressed.connect(_mode.spawn_free.bind("proton"))
	%SpawnNBtn.pressed.connect(_mode.spawn_free.bind("neutron"))
	%SpawnEBtn.pressed.connect(_mode.spawn_free.bind("electron"))

	%SubstepsSlider.value_changed.connect(_on_substeps_changed)
	%SubstepsSpinBox.value_changed.connect(_on_substeps_changed)

func _on_substeps_changed(value: float) -> void:
	_mode.substeps_per_frame = int(value)
	%SubstepsSlider.set_value_no_signal(value)
	%SubstepsSpinBox.set_value_no_signal(value)

func _process(_delta: float) -> void:
	if _mode == null or _mode.atom_sim == null:
		return
	var sim: Node = _mode.atom_sim

	var paused: bool = _mode.paused
	%StatusLabel.text = "paused" if paused else "running"
	%PauseBtn.text = "Resume" if paused else "Pause"
	%PauseLabel.text = "PAUSED" if paused else ""
	%TimeLabel.text = "%.3f" % float(sim.get_time())
	%ParticlesLabel.text = "%d" % int(sim.get_particle_count())
	%KELabel.text = String.num_scientific(float(sim.get_total_kinetic_energy()))

	var oi: PackedFloat32Array = sim.get_orbit_info()
	if oi.size() >= 6 and oi[0] > 0.0:
		%OrbitRLabel.text = "%.3f  /  %.1f°" % [oi[0], oi[5]]
		%OrbitVLabel.text = "%.3f  /  %+.3f" % [oi[1], oi[2]]
	else:
		%OrbitRLabel.text = "—"
		%OrbitVLabel.text = "—"

	var fb: PackedFloat32Array = sim.get_force_breakdown()
	if fb.size() >= 8:
		%FGravLabel.text = "%s / %s" % [
			String.num_scientific(fb[0]), String.num_scientific(fb[1])]
		%FIntakeLabel.text = "%s / %s" % [
			String.num_scientific(fb[2]), String.num_scientific(fb[3])]
	else:
		%FGravLabel.text = "—"
		%FIntakeLabel.text = "—"

	if int(sim.get_particle_count()) >= 2:
		%PairLabel.text = "%.4f" % float(sim.get_pair_distance(0, 1))
	else:
		%PairLabel.text = "—"

	%FpsLabel.text = "%d" % Engine.get_frames_per_second()

	# Mirror toggles that can also change via hotkeys.
	%CloudsBtn.set_pressed_no_signal(_mode.show_clouds)
	if int(%SubstepsSlider.value) != _mode.substeps_per_frame:
		%SubstepsSlider.set_value_no_signal(_mode.substeps_per_frame)
		%SubstepsSpinBox.set_value_no_signal(_mode.substeps_per_frame)
