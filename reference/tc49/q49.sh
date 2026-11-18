#!/bin/bash
# Doc 15.49 (ESTELLE half): near-exit collision count for the neutron's beam, p-n frozen at d 3 and
# 2.05, local pairing off and on, seed 49344 local pairing off; ISABELLE runs local pairing on, seed 777 (q49_*_lp1_s777). First pass printed no beam lines (len check bug), kept in old_noprint/
export PATH=/usr/bin:/bin:$PATH
cd /c/Users/nnt04/Project/charge-field-modeler/rust
B=target/fast/deps/charge_field_modeler-8ba177252dd33c00.exe
O=../reference/tc49
P=$O/progress.txt
echo "START $(date)" > $P
export CFM_D_INTAKE=ball CFM_D_COPIES=64 CFM_D_POLAR=through CFM_D_BEAM=0.222 CFM_C_LAMBDA=0.3 CFM_C_H=0.25 CFM_C_DT=0.06 CFM_F_RD=6
export CFM_D_FREEZE=1 CFM_D_FIXT=150 CFM_D_CALWARMT=150 CFM_D_CALT=100 CFM_D_STEPS=1667 CFM_D_LOG=167
export CFM_D_FSELF="proton:-106.9,-12.0,22.6,1.5,1.8,-4.2;neutron:1.5,5.8,-5.0,-3.3,1.5,-5.8"
export CFM_F_PAIR=pn
run() { $B --ignored --nocapture --exact gpu::dsmc::shell_gates::report_deuteron > $O/$1.txt 2>&1; echo "$1 done $(date)" >> $P; }
for d in 3 2.05; do for lp in 0; do CFM_D_SEP0=$d CFM_D_LOCALPAIR=$lp run pn_${d}_lp${lp}_s49344; done; done
echo "ALL DONE $(date)" >> $P
