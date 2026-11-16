#!/bin/bash
# Doc 15.44: how robust is the on-axis He3 (p-n-p, 0.685 neutron, through-release, channels)?
# 1-2: frozen at centres 2.3 and 2.6 (where does the force change sign between 2.05 and 3?)
# 3:   frozen at 2.05 with the +z proton's axis flipped (net-force check; antiparallel ends)
# 4:   free run released from 2.3 (does it fall back to contact or come apart?)
cd /c/Users/nnt04/Project/charge-field-modeler/rust
BIN=/c/Users/nnt04/Project/charge-field-modeler/rust/target/fast/deps/charge_field_modeler-8ba177252dd33c00.exe
O=/c/Users/nnt04/Project/charge-field-modeler/reference/ch47
export CFM_F_RD=6 CFM_D_LOG=500 CFM_D_INTAKE_N=deflect0.915 CFM_D_POLAR=through CFM_D_HOLE=0.3333
export CFM_F_PAIR=pnp CFM_D_LAYOUT=custom
export CFM_D_FSELF="proton:2.6,-4.0,45.4,0,0,0;neutron:4.1,6.4,29.1,-9.0,-0.7,7.0"
run() { "$BIN" --ignored --nocapture --exact gpu::dsmc::shell_gates::report_deuteron > "$O/$1.txt" 2>&1; echo "$1 done $(date)" >> "$O/progress.txt"; }
CFM_D_FREEZE=1 CFM_D_STEPS=7500 CFM_D_POS="0,0,-2.3;0,0,0;0,0,2.3" run he3_frozen_230
CFM_D_FREEZE=1 CFM_D_STEPS=7500 CFM_D_POS="0,0,-2.6;0,0,0;0,0,2.6" run he3_frozen_260
CFM_D_FREEZE=1 CFM_D_STEPS=7500 CFM_D_POS="0,0,-2.05;0,0,0;0,0,2.05" CFM_D_AXES="0,0,1;0,0,1;0,0,-1" run he3_frozen_flip
CFM_D_STEPS=30000 CFM_D_MASS=7e6 CFM_D_INERTIA=2.8e6 CFM_D_LSPIN=8.75e5 CFM_D_POS="0,0,-2.3;0,0,0;0,0,2.3" run he3_free_from_230
echo "ALL DONE $(date)" >> "$O/progress.txt"
