#!/bin/bash
# Doc 15.44: two alpha ideas with the current neutron (0.685, parallel through-release) and
# proton channels 1/3. Compact alpha: protons z +-1.9 (same spin), neutrons x +-1.0.
# 1: Mathis's sandwich (neutron axes +-z): the beams (radius 1/3 at x +-1) straddle the proton rims.
# 2: neutrons on their side, in line along x, head to tail (pole to pole, beams leave along +-x).
# 3: free run of 2 (rotation on, 600 tu).
# Waits for the ch47 queue to finish first.
cd /c/Users/nnt04/Project/charge-field-modeler/rust
until grep -q "ALL DONE" ../reference/ch47/progress.txt 2>/dev/null; do sleep 60; done
BIN=/c/Users/nnt04/Project/charge-field-modeler/rust/target/fast/deps/charge_field_modeler-8ba177252dd33c00.exe
O=/c/Users/nnt04/Project/charge-field-modeler/reference/ch48
export CFM_F_RD=6 CFM_D_LOG=500 CFM_D_INTAKE_N=deflect0.915 CFM_D_POLAR=through CFM_D_HOLE=0.3333
export CFM_F_PAIR=pnnp CFM_D_LAYOUT=custom
export CFM_D_FSELF="proton:2.6,-4.0,45.4,0,0,0;neutron:4.1,6.4,29.1,-9.0,-0.7,7.0"
run() { "$BIN" --ignored --nocapture --exact gpu::dsmc::shell_gates::report_deuteron > "$O/$1.txt" 2>&1; echo "$1 done $(date)" >> "$O/progress.txt"; }
P="0,0,-1.9;-1,0,0;1,0,0;0,0,1.9"
CFM_D_FREEZE=1 CFM_D_STEPS=7500 CFM_D_POS="$P" CFM_D_AXES="0,0,1;0,0,1;0,0,-1;0,0,1" run sandwich_through
CFM_D_FREEZE=1 CFM_D_STEPS=7500 CFM_D_POS="$P" CFM_D_AXES="0,0,1;1,0,0;1,0,0;0,0,1" run sideways_inline
CFM_D_STEPS=30000 CFM_D_MASS=7e6 CFM_D_INERTIA=2.8e6 CFM_D_LSPIN=8.75e5 CFM_D_POS="0,0,-1.9;-1.02,0,0;1.02,0,0;0,0,1.9" CFM_D_AXES="0,0,1;1,0,0;1,0,0;0,0,1" run sideways_inline_free
echo "ALL DONE $(date)" >> "$O/progress.txt"
