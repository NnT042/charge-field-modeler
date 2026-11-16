#!/bin/bash
# Doc 15.44 free runs: on-axis p-n-p and p-n-n-p from contact (centres 2.05), 0.685 neutron,
# parallel through-release, proton channel 1/3. Free in 3D (custom layout), rotation on
# (settings of the 15.34/15.36 free runs): 250 tu held, then 600 tu free. Two seeds each.
cd /c/Users/nnt04/Project/charge-field-modeler/rust
BIN=/c/Users/nnt04/Project/charge-field-modeler/rust/target/fast/deps/charge_field_modeler-8ba177252dd33c00.exe
O=/c/Users/nnt04/Project/charge-field-modeler/reference/ch46
export CFM_F_RD=6 CFM_D_STEPS=30000 CFM_D_LOG=500 CFM_D_INTAKE_N=deflect0.915 CFM_D_POLAR=through CFM_D_HOLE=0.3333
export CFM_D_MASS=7e6 CFM_D_INERTIA=2.8e6 CFM_D_LSPIN=8.75e5 CFM_D_LAYOUT=custom
export CFM_D_FSELF="proton:2.6,-4.0,45.4,0,0,0;neutron:4.1,6.4,29.1,-9.0,-0.7,7.0"
run() { "$BIN" --ignored --nocapture --exact gpu::dsmc::shell_gates::report_deuteron > "$O/$1.txt" 2>&1; echo "$1 done $(date)" >> "$O/progress.txt"; }
HE3="0,0,-2.05;0,0,0;0,0,2.05"
ALPHA="0,0,-3.075;0,0,-1.025;0,0,1.025;0,0,3.075"
CFM_C_SEED=49344 CFM_F_PAIR=pnp  CFM_D_POS="$HE3"   run he3_free_s49344
CFM_C_SEED=49344 CFM_F_PAIR=pnnp CFM_D_POS="$ALPHA" run pnnp_free_s49344
CFM_C_SEED=777   CFM_F_PAIR=pnp  CFM_D_POS="$HE3"   run he3_free_s777
CFM_C_SEED=777   CFM_F_PAIR=pnnp CFM_D_POS="$ALPHA" run pnnp_free_s777
echo "ALL DONE $(date)" >> "$O/progress.txt"
