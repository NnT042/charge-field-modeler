#!/bin/bash
# Doc 15.44: 0.685 neutron (deflect 0.915), parallel through-release beam (radius 1/3),
# aimed down the proton's through-channel (hole 1/3). Frozen screens, 15.41 settings.
cd /c/Users/nnt04/Project/charge-field-modeler/rust
BIN=/c/Users/nnt04/Project/charge-field-modeler/rust/target/fast/deps/charge_field_modeler-8ba177252dd33c00.exe
O="/c/Users/nnt04/AppData/Local/Temp/claude/63a7c279081ac85c/b395b8e1-54e5-4c45-a8ec-a703a023d0dd/scratchpad/ch44"
export CFM_F_RD=6 CFM_D_FREEZE=1 CFM_D_STEPS=7500 CFM_D_LOG=500 CFM_D_INTAKE_N=deflect0.915 CFM_D_POLAR=through
PH="proton:2.6,-4.0,45.4,0,0,0"
PN="proton:-10.4,-14.1,63.8,9.5,14.0,-7.6"
run() { "$BIN" --ignored --nocapture --exact gpu::dsmc::shell_gates::report_deuteron > "$O/$1.txt" 2>&1; echo "$1 done $(date)" >> "$O/progress.txt"; }
# 1: on-axis deuteron, centres 3, proton WITH channel (calibrates the lone neutron)
CFM_F_PAIR=pn CFM_D_SEP0=3 CFM_D_HOLE=0.3333 CFM_D_FSELF="$PH" run pn_axis_hole
F=$(grep "^self-force neutron" "$O/pn_axis_hole.txt" | head -1 | sed -E 's/.*\(([^)]*)\).*/\1/; s/[ +]//g')
T=$(grep "^self-torque neutron" "$O/pn_axis_hole.txt" | head -1 | sed -E 's/.*\(([^)]*)\).*/\1/; s/[ +]//g')
NS="neutron:$F,$T"
# 2: same, proton WITHOUT channel (control: the beam is absorbed)
CFM_F_PAIR=pn CFM_D_SEP0=3 CFM_D_FSELF="$PN;$NS" run pn_axis_nohole
# 3: He3 hub p-n-p on axis, protons with channel
CFM_F_PAIR=pnp CFM_D_SEPS=2.05,2.05 CFM_D_HOLE=0.3333 CFM_D_FSELF="$PH;$NS" run he3_axis_hole
echo "ALL DONE $(date)" >> "$O/progress.txt"
