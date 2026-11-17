#!/bin/bash
# Doc 15.45: trace bodies (64 copies, surface emission, lambda 0.3) with the neutron's trapped
# photons leaving by the pole they were heading toward (CFM_D_POLAR=through, user's §15.32 rule).
cd /c/Users/nnt04/Project/charge-field-modeler/rust
B=target/fast/deps/charge_field_modeler-8ba177252dd33c00.exe
O=../reference/tc45
export CFM_D_INTAKE=ball CFM_D_COPIES=64 CFM_D_POLAR=through CFM_C_LAMBDA=0.3 CFM_C_H=0.25 CFM_C_DT=0.06 CFM_F_RD=6
export CFM_D_FREEZE=1 CFM_D_FIXT=150 CFM_D_CALWARMT=150 CFM_D_CALT=100 CFM_D_STEPS=1667 CFM_D_LOG=167
run() { $B --ignored --nocapture --exact gpu::dsmc::shell_gates::report_deuteron > $O/$1.txt 2>&1; echo "$1 done $(date)" >> $O/progress_through.txt; }
CFM_F_PAIR=pn CFM_D_SEP0=3 run tr_pn_axis_3_through
F=$(grep "^self-force neutron" $O/tr_pn_axis_3_through.txt | head -1 | sed -E 's/.*\(([^)]*)\).*/\1/; s/[ +]//g'); T=$(grep "^self-torque neutron" $O/tr_pn_axis_3_through.txt | head -1 | sed -E 's/.*\(([^)]*)\).*/\1/; s/[ +]//g')
FS="proton:-106.9,-12.0,22.6,1.5,1.8,-4.2;neutron:$F,$T"
echo "FSELF $FS" >> $O/progress_through.txt
CFM_F_PAIR=pnp CFM_D_SEPS=2.05,2.05 CFM_D_FSELF="$FS" run tr_he3_axis_205_through
echo "ALL DONE $(date)" >> $O/progress_through.txt
