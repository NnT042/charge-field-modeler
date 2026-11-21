#!/bin/bash
# First rented queue (Doc 15.57). Runs JOBS at a time on the one GPU; each job is a tc56-style
# report_deuteron run. Start inside tmux:  tmux new -d -s q 'bash q57.sh'
# Results: ~/out/*.txt, progress in ~/out/progress.txt.
JOBS=${JOBS:-4}
B=$HOME/cfm_test
O=$HOME/out; mkdir -p $O; P=$O/progress.txt
cd $HOME/cfm/rust
echo "START $(date) jobs $JOBS" >> $P
export CFM_D_INTAKE=ball CFM_D_COPIES=64 CFM_D_POLAR=through CFM_D_BEAM=0.222 CFM_F_RD=6 CFM_C_SEED=49344
export CFM_D_FIXT=150 CFM_D_CALWARMT=150 CFM_D_CALT=100 CFM_D_LOG=167 CFM_C_DT=0.06 CFM_C_LAMBDA=0.3 CFM_C_H=0.25
LANE="CFM_D_PIPE=pn CFM_D_PIPE_DEG=90 CFM_D_PIPE_LANE=1"
ROT="CFM_D_INERTIA=2.8e6 CFM_D_LSPIN=8.75e5"
TRI="CFM_D_LAYOUT=custom CFM_D_POS=0,0,-2.05;0,0,0;0,0,2.05 CFM_D_MASS=7e6 CFM_D_STEPS=10000"
HI="CFM_C_NA=16000 CFM_C_H=0.125"
run() { n=$1; shift; ( export "$@"; $B --ignored --nocapture --exact gpu::dsmc::shell_gates::report_deuteron > $O/$n.txt 2>&1; echo "$n done $(date) exit $?" >> $P ); }
q() { while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do wait -n; done; run "$@" & }
# 1. Free n-p-n with the lane (strongest frozen pull, +400): three seeds.
q npn_free_lane_s49344 $LANE $ROT $TRI CFM_F_PAIR=npn
q npn_free_lane_s1     $LANE $ROT $TRI CFM_F_PAIR=npn CFM_C_SEED=1
q npn_free_lane_s2     $LANE $ROT $TRI CFM_F_PAIR=npn CFM_C_SEED=2
# 2. Free p-n-p with the lane: two more seeds (tc56 came apart at ~100 tu on 49344).
q pnp_free_lane_s1     $LANE $ROT $TRI CFM_F_PAIR=pnp CFM_C_SEED=1
q pnp_free_lane_s2     $LANE $ROT $TRI CFM_F_PAIR=pnp CFM_C_SEED=2
# 3. Windsock without the gyroscope (L_spin 0): does the 10-degree tilt swing back when the torque is not turned into precession?
q pn_tilt10_lane_nogyro $LANE CFM_D_INERTIA=2.8e6 CFM_D_LSPIN=0 CFM_F_PAIR=pn CFM_D_LAYOUT=custom "CFM_D_POS=0,0,-1.025;0,0,1.025" "CFM_D_AXES=0.1736,0,0.9848;0,0,1" CFM_D_MASS=1e15 CFM_D_STEPS=5000
# 4. Hi-res free runs (lane +219 at h 0.125): does the wider/deeper well hold the triple?
q npn_free_lane_hires  $LANE $ROT $TRI $HI CFM_F_PAIR=npn
q pnp_free_lane_hires  $LANE $ROT $TRI $HI CFM_F_PAIR=pnp
wait
echo "ALL DONE $(date)" >> $P
