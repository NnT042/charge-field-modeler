#!/bin/bash
# Doc 15.61 (ESTELLE half): headwind, lone proton / antiproton moving through the masspath photon mix
# (hit box + mass, even mix, R_D 6, pre-filled). ISABELLE runs q61_a..f (gates, p at V 0.01-0.05, p tilt 30 at V 0.02).
export PATH=/usr/bin:/bin:$PATH
cd /c/Users/nnt04/Project/charge-field-modeler/rust
B=/c/Users/nnt04/AppData/Local/Temp/claude/63a7c279081ac85c/b395b8e1-54e5-4c45-a8ec-a703a023d0dd/scratchpad/cfm61.exe
O=../reference/tc61
P=$O/progress_estelle.txt
echo "START $(date)" > $P
export CFM_C_SEED=49344 CFM_D_INTAKE=ball CFM_D_COPIES=64 CFM_D_POLAR=through CFM_D_BEAM=0.222 CFM_C_DT=0.06 CFM_F_RD=6
export CFM_D_CALWARMT=150 CFM_D_CALT=100 CFM_C_LAMBDA=0.3 CFM_C_H=0.25
export CFM_P_MIX=1,1,1,1 CFM_P_INNER=random CFM_P_SIZE=1 CFM_P_MASS=1 CFM_P_PREFILL=1 CFM_P_NORM=masspath CFM_P_CELLS=1 CFM_D_LONEONLY=1
go() { n=$1; shift; ( export "$@"; $B --ignored --nocapture --exact gpu::dsmc::shell_gates::report_deuteron > $O/$n.txt 2>&1; echo "$n done $(date)" >> $P ); }
go q61_g_p_vneg020 CFM_F_PAIR=p CFM_D_WIND=0,0,0.02
go q61_h_p_tilt90_v020 CFM_F_PAIR=p CFM_D_WIND=0,0,-0.02 CFM_D_LONEAXIS=1,0,0
go q61_i_a_v020 CFM_F_PAIR=a CFM_D_WIND=0,0,0.02
go q61_j_p_tilt30_v010 CFM_F_PAIR=p CFM_D_WIND=0,0,-0.01 CFM_D_LONEAXIS=0.5,0,0.8660254
echo "ALL DONE $(date)" >> $P
