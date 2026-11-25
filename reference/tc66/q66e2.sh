#!/bin/bash
# Doc 15.66 (ESTELLE, round 2): rigid maze spin, seed 2 and seed 49344 for 600 tu in 6 blocks.
# ISABELLE runs q66_a..c (gates; rigid maze spin seed 49344 and seed 1).
export PATH=/usr/bin:/bin:$PATH
cd /c/Users/nnt04/Project/charge-field-modeler/rust
B=/c/Users/nnt04/AppData/Local/Temp/claude/63a7c279081ac85c/b395b8e1-54e5-4c45-a8ec-a703a023d0dd/scratchpad/cfm66b.exe
O=../reference/tc66
P=$O/progress_estelle2.txt
echo "START $(date)" > $P
export CFM_C_SEED=49344 CFM_D_INTAKE=ball CFM_D_COPIES=64 CFM_D_POLAR=through CFM_D_BEAM=0.222 CFM_C_DT=0.06 CFM_F_RD=6
export CFM_D_CALWARMT=150 CFM_D_CALT=300 CFM_C_LAMBDA=0.3 CFM_C_H=0.25
export CFM_P_MIX=1,1,1,1 CFM_P_INNER=random CFM_P_SIZE=1 CFM_P_MASS=1 CFM_P_PREFILL=1 CFM_P_NORM=masspath CFM_D_LONEONLY=1 CFM_F_PAIR=p
export CFM_D_CALBLOCKS=3 CFM_D_MAZEPHASE=random
go() { n=$1; shift; ( export "$@"; $B --ignored --nocapture --exact gpu::dsmc::shell_gates::report_deuteron > $O/$n.txt 2>&1; echo "$n done $(date)" >> $P ); }
go q66_h_rigid_seed2 CFM_C_SEED=2 CFM_D_MAZESPIN=rigid
go q66_i_rigid_base_long CFM_D_MAZESPIN=rigid CFM_D_CALT=600 CFM_D_CALBLOCKS=6
echo "ALL DONE $(date)" >> $P
