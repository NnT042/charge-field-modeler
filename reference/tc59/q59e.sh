#!/bin/bash
# Doc 15.59 (ESTELLE): photon mix at bodies. Frozen p-n at 2.05 (also prints the lone proton and lone
# neutron), then the free H atom (p-e). Standard tc56 config. Even born mix (CFM_P_MIX 1,1,1,1).
# ISABELLE (C:\cfm\queue, q59_*) runs the hit-box versions (CFM_P_SIZE=1), which are ~30x slower.
export PATH=/usr/bin:/bin:$PATH
cd /c/Users/nnt04/Project/charge-field-modeler/rust
B=/c/Users/nnt04/AppData/Local/Temp/claude/63a7c279081ac85c/b395b8e1-54e5-4c45-a8ec-a703a023d0dd/scratchpad/cfm59.exe
O=../reference/tc59
P=$O/progress_estelle.txt
echo "START $(date)" > $P
export CFM_D_INTAKE=ball CFM_D_COPIES=64 CFM_D_POLAR=through CFM_D_BEAM=0.222 CFM_F_RD=6 CFM_C_SEED=49344
export CFM_D_FIXT=150 CFM_D_CALWARMT=150 CFM_D_CALT=100 CFM_D_LOG=167 CFM_C_DT=0.06 CFM_C_LAMBDA=0.3 CFM_C_H=0.25
go() { n=$1; shift; ( export "$@"; run $n ); }
run() { $B --ignored --nocapture --exact gpu::dsmc::shell_gates::report_deuteron > $O/$1.txt 2>&1; echo "$1 done $(date)" >> $P; }
PN="CFM_F_PAIR=pn CFM_D_FREEZE=1 CFM_D_STEPS=1667 CFM_D_SEP0=2.05"
EQ="CFM_P_MIX=1,1,1,1"
go pn_k0 $PN
for inner in plus minus random; do go pn_eq_$inner $PN $EQ CFM_P_INNER=$inner; done
for inner in plus minus random; do go pn_eq_${inner}_mass $PN $EQ CFM_P_INNER=$inner CFM_P_MASS=1; done
PE="CFM_F_PAIR=pe CFM_D_LAYOUT=custom CFM_D_POS=0,0,0;0,0,2.5 CFM_D_MASS=7e6 CFM_D_STEPS=10000"
go pe_k0 $PE
go pe_eq_random_mass $PE $EQ CFM_P_INNER=random CFM_P_MASS=1
echo "ALL DONE $(date)" >> $P
