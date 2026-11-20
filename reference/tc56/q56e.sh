#!/bin/bash
# Doc 15.56 overnight (ESTELLE): hi-res pilot, free He3 with the lane, windsock tilt, bubble maps.
# Trace bodies as tc55 (ball intake, 64 copies, through release, beam 0.222, lambda 0.3, h 0.25, R_D 6, seed 49344).
# Every run calibrates its own self-forces. ISABELLE runs the 3-body screens and the electron recheck.
export PATH=/usr/bin:/bin:$PATH
cd /c/Users/nnt04/Project/charge-field-modeler/rust
B=/c/Users/nnt04/AppData/Local/Temp/claude/63a7c279081ac85c/b395b8e1-54e5-4c45-a8ec-a703a023d0dd/scratchpad/cfm56.exe
O=../reference/tc56
P=$O/progress_estelle.txt
echo "START $(date)" > $P
export CFM_D_INTAKE=ball CFM_D_COPIES=64 CFM_D_POLAR=through CFM_D_BEAM=0.222 CFM_F_RD=6 CFM_F_PAIR=pn CFM_C_SEED=49344
export CFM_D_FIXT=150 CFM_D_CALWARMT=150 CFM_D_CALT=100 CFM_D_LOG=167 CFM_C_DT=0.06 CFM_C_LAMBDA=0.3 CFM_C_H=0.25
go() { n=$1; shift; ( export "$@"; run $n ); }
run() { $B --ignored --nocapture --exact gpu::dsmc::shell_gates::report_deuteron > $O/$1.txt 2>&1; echo "$1 done $(date)" >> $P; }
LANE="CFM_D_PIPE=pn CFM_D_PIPE_DEG=90 CFM_D_PIPE_LANE=1"
ROT="CFM_D_INERTIA=2.8e6 CFM_D_LSPIN=8.75e5"
# 1. Hi-res pilot: twice the photons (n_a 16000, ~14.5M), cells half as wide (h 0.125), same mean free path. Frozen p-n, lane.
go hires_pn_2.05_lane $LANE CFM_C_NA=16000 CFM_C_H=0.125 CFM_D_FREEZE=1 CFM_D_STEPS=1667 CFM_D_SEP0=2.05
# 2. Free He3 (p-n-p) from contact, 3D, rotation on, lane on both classes, 600 tu (ch46 ran this with the old sponge bodies: held).
go he3_free_lane $LANE $ROT CFM_F_PAIR=pnp CFM_D_LAYOUT=custom CFM_D_POS="0,0,-2.05;0,0,0;0,0,2.05" CFM_D_MASS=7e6 CFM_D_STEPS=10000
# 3. Windsock: p-n at 2.05 with the lane, centres held (huge mass), proton axis tilted 10 deg; does the axis swing back?
go pn_tilt10_lane $LANE $ROT CFM_D_LAYOUT=custom CFM_D_POS="0,0,-1.025;0,0,1.025" CFM_D_AXES="0.1736,0,0.9848;0,0,1" CFM_D_MASS=1e15 CFM_D_STEPS=5000
# 4. Bubble maps (frozen): merged bubble at contact with and without the lane, and a pair 4 apart.
go map_pn_2.05_lane $LANE CFM_D_FREEZE=1 CFM_D_STEPS=1667 CFM_D_SEP0=2.05 CFM_D_MAP=$O/map_pn_2.05_lane.json
go map_pn_2.05 CFM_D_FREEZE=1 CFM_D_STEPS=1667 CFM_D_SEP0=2.05 CFM_D_MAP=$O/map_pn_2.05.json
go map_pn_4 CFM_D_FREEZE=1 CFM_D_STEPS=1667 CFM_D_SEP0=4 CFM_D_MAP=$O/map_pn_4.json
# 5. Hi-res pilot of the best channel rule (edge-on, eta 0.3).
go hires_pn_2.05_edge03 CFM_D_PIPE=pn CFM_D_PIPE_DEG=90 CFM_D_PIPE_SPIN=1 CFM_D_PIPE_ETA=0.3 CFM_D_PIPE_EDGE=1 CFM_C_NA=16000 CFM_C_H=0.125 CFM_D_FREEZE=1 CFM_D_STEPS=1667 CFM_D_SEP0=2.05
echo "ALL DONE $(date)" >> $P
