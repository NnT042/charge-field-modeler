#!/bin/bash
# Doc 15.45 pair screens: trace bodies (Intake::Ball, maze, surface emission), 64 trace copies,
# frozen, lambda 0.3, n_a 8000, h 0.25, dt 0.06, R_D 6; 150 tu settle + 100 tu measured.
cd /c/Users/nnt04/Project/charge-field-modeler/rust
B=target/fast/deps/charge_field_modeler-8ba177252dd33c00.exe
O=../reference/tc45
export CFM_D_INTAKE=ball CFM_D_COPIES=64 CFM_C_LAMBDA=0.3 CFM_C_H=0.25 CFM_C_DT=0.06 CFM_F_RD=6
export CFM_D_FREEZE=1 CFM_D_FIXT=150 CFM_D_CALWARMT=150 CFM_D_CALT=100 CFM_D_STEPS=1667 CFM_D_LOG=167
run() { $B --ignored --nocapture --exact gpu::dsmc::shell_gates::report_deuteron > $O/$1.txt 2>&1; echo "$1 done $(date)" >> $O/progress_pairs.txt; }
sf() { F=$(grep "^self-force $1" $O/$2.txt | head -1 | sed -E 's/.*\(([^)]*)\).*/\1/; s/[ +]//g'); T=$(grep "^self-torque $1" $O/$2.txt | head -1 | sed -E 's/.*\(([^)]*)\).*/\1/; s/[ +]//g'); echo "$1:$F,$T"; }
# 1: deuteron pole-on at 3 (calibrates both classes)
CFM_F_PAIR=pn CFM_D_SEP0=3 run tr_pn_axis_3
FS="$(sf proton tr_pn_axis_3);$(sf neutron tr_pn_axis_3)"
echo "FSELF $FS" >> $O/progress_pairs.txt
# 2: proton pair pole-on at 3
CFM_F_PAIR=pp CFM_D_SEP0=3 CFM_D_FSELF="$FS" run tr_pp_axis_3
# 3: deuteron side by side at 3
CFM_F_PAIR=pn CFM_F_GEOM=side CFM_D_SEP0=3 CFM_D_FSELF="$FS" run tr_pn_side_3
# 4: He3 on axis at contact
CFM_F_PAIR=pnp CFM_D_SEPS=2.05,2.05 CFM_D_FSELF="$FS" run tr_he3_axis_205
echo "ALL DONE $(date)" >> $O/progress_pairs.txt
