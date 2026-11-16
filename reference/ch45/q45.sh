#!/bin/bash
# Doc 15.44 follow-up: same settings as reference/ch44/q44.sh (0.685 neutron, parallel
# through-release, proton channel 1/3, frozen). Separates distance from layout, and tries
# two on-axis alphas.
cd /c/Users/nnt04/Project/charge-field-modeler/rust
BIN=/c/Users/nnt04/Project/charge-field-modeler/rust/target/fast/deps/charge_field_modeler-8ba177252dd33c00.exe
O=/c/Users/nnt04/Project/charge-field-modeler/reference/ch45
export CFM_F_RD=6 CFM_D_FREEZE=1 CFM_D_STEPS=7500 CFM_D_LOG=500 CFM_D_INTAKE_N=deflect0.915 CFM_D_POLAR=through CFM_D_HOLE=0.3333
export CFM_D_FSELF="proton:2.6,-4.0,45.4,0,0,0;neutron:4.1,6.4,29.1,-9.0,-0.7,7.0"
run() { "$BIN" --ignored --nocapture --exact gpu::dsmc::shell_gates::report_deuteron > "$O/$1.txt" 2>&1; echo "$1 done $(date)" >> "$O/progress.txt"; }
# 1: deuteron at contact distance (the ch44 deuteron was at 3, the He3 at 2.05)
CFM_F_PAIR=pn CFM_D_SEP0=2.05 run pn_axis_hole_205
# 2: on-axis alpha, neutrons in the middle
CFM_F_PAIR=pnnp CFM_D_SEPS=2.05,2.05,2.05 run alpha_pnnp
# 3: on-axis alpha, protons in the middle (hole to hole), neutrons outside
CFM_F_PAIR=nppn CFM_D_SEPS=2.05,2.05,2.05 run alpha_nppn
# 4: He3 at centres 3 (distance check the other way)
CFM_F_PAIR=pnp CFM_D_SEPS=3,3 run he3_axis_hole_3
echo "ALL DONE $(date)" >> "$O/progress.txt"
