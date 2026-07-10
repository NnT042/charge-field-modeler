extends SceneTree
#
# TEMP session-31 probe: reproduce the user-reported carbon dissolution in
# RigidAlpha using the REAL app scene (HUD + tuning panel + VFX included).
# Run: godot --headless --path godot -s res://scripts/dev_probe_dissolve.gd
#

var frames := 0
var mode: Node = null
var d0 := {}

# Particle-index pairs for the BARE 3-STACK ("tri_alpha" since session 32
# — the structure the session-31 dissolution reports were filed against;
# carbon's real shape is now 2 plugged alphas, see preset_alphas).
const PAIRS := [[0, 4], [4, 8], [0, 8], [0, 1], [8, 9]]

func _initialize() -> void:
	var packed: PackedScene = load("res://scenes/atom_mode.tscn")
	mode = packed.instantiate()
	root.add_child(mode)

func _sample() -> Dictionary:
	var d := {}
	for p in PAIRS:
		d["%d-%d" % [p[0], p[1]]] = float(mode.atom_sim.get_pair_distance(p[0], p[1]))
	return d

func _process(_delta: float) -> bool:
	frames += 1
	if frames == 5:
		mode.spawn_scenario("tri_alpha")
	elif frames == 60:
		mode.set_nucleus_dynamics(1)
		d0 = _sample()
		print("[probe] toggled RigidAlpha, d0=", d0)
	elif frames > 60 and frames % 120 == 0:
		var d := _sample()
		var worst := 0.0
		for k in d:
			var rel: float = abs(d[k] - d0[k]) / max(d0[k], 1e-9)
			worst = max(worst, rel)
		print("[probe] frame=%d worst_drift=%.1f%% d=%s" % [frames, worst * 100.0, str(d)])
		if worst > 3.0:
			print("[probe] DISSOLVED")
			quit(0)
			return true
	# 10000 frames x 100 substeps = 1M steps — past the deepest pre-fix
	# whirl-instability horizon (825k, seed k=7; session 32).
	if frames >= 10000:
		print("[probe] done, stable")
		quit(0)
		return true
	return false
