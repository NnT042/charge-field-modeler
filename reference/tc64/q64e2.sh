#!/bin/bash
# Doc 15.64 (ESTELLE half): when does the lopsided state appear, and does the polar stream hold it.
# ISABELLE runs q64_j..l (emergence base; emergence roll 180; polar axial).
export PATH=/usr/bin:/bin:$PATH
cd /c/Users/nnt04/Project/charge-field-modeler/rust
B=/c/Users/nnt04/AppData/Local/Temp/claude/63a7c279081ac85c/b395b8e1-54e5-4c45-a8ec-a703a023d0dd/scratchpad/cfm63.exe
O=../reference/tc64
P=$O/progress_estelle.txt
echo "START2 $(date)" >> $P
export CFM_C_SEED=49344 CFM_D_INTAKE=ball CFM_D_COPIES=64 CFM_D_POLAR=through CFM_D_BEAM=0.222 CFM_C_DT=0.06 CFM_F_RD=6
export CFM_D_CALWARMT=150 CFM_D_CALT=300 CFM_C_LAMBDA=0.3 CFM_C_H=0.25
export CFM_P_MIX=1,1,1,1 CFM_P_INNER=random CFM_P_SIZE=1 CFM_P_MASS=1 CFM_P_PREFILL=1 CFM_P_NORM=masspath CFM_D_LONEONLY=1 CFM_F_PAIR=p
export CFM_D_CALBLOCKS=3
go() { n=$1; shift; ( export "$@"; $B --ignored --nocapture --exact gpu::dsmc::shell_gates::report_deuteron > $O/$n.txt 2>&1; echo "$n done $(date)" >> $P ); }
go q64_p_mazeiso_roll180 CFM_D_MAZE=iso CFM_D_LONEROLL=180
echo "ALL DONE2 $(date)" >> $P
