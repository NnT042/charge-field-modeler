extends SceneTree
#
# Headless integration smoke for atom mode that actually exercises the
# scenario paths (presets -> rigid groups -> composite skins, h2 -> bond
# bridges), which the plain scene smoke never leaves "protons" for.
#
# Run:
#   godot --headless --path godot -s res://scripts/dev_smoke_atom.gd
#
# Pass = exit 0, prints "[smoke] OK", no SCRIPT ERROR lines.

var frames := 0
var mode: Node = null
var failed := false

func _initialize() -> void:
	var packed: PackedScene = load("res://scenes/atom_mode.tscn")
	mode = packed.instantiate()
	root.add_child(mode)

func _process(_delta: float) -> bool:
	frames += 1
	match frames:
		5:
			mode.spawn_scenario("alpha")
		10:
			_check_groups(1, "alpha")
			mode.spawn_scenario("h2_bond")
		15:
			var bonds: PackedInt32Array = mode.atom_sim.get_bond_pairs()
			if bonds.size() < 2:
				push_error("[smoke] h2_bond: no molecular bond detected")
				failed = true
		20:
			mode.spawn_scenario("oxygen")
		25:
			_check_groups(1, "oxygen")
			mode.spawn_scenario("protons")
		30:
			_check_groups(0, "protons")
	if frames >= 35:
		if failed:
			print("[smoke] FAILED")
			quit(1)
		else:
			print("[smoke] OK")
			quit(0)
		return true
	return false

func _check_groups(expected: int, label: String) -> void:
	var gcount: int = mode.atom_sim.get_group_count()
	var skins: int = mode.group_skins.size()
	if gcount != expected or skins != expected:
		push_error("[smoke] %s: groups=%d skins=%d expected=%d" % [label, gcount, skins, expected])
		failed = true
