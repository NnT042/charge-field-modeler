#!/bin/bash
# Hit-box resolution ladder (docs/VAST_HITBOX_PLAN.md, doc 15.59). Runs JOBS at a time on the one GPU.
# Start inside tmux:  tmux new -d -s q 'RUNG=1 bash q59.sh'   (RUNG = 1, 2 or 3)
# Results: ~/out/*.txt, progress in ~/out/progress.txt.
JOBS=${JOBS:-4}
RUNG=${RUNG:-1}
B=$HOME/cfm_test
O=$HOME/out; mkdir -p $O; P=$O/progress.txt
cd $HOME/cfm/rust
echo "START $(date) rung $RUNG jobs $JOBS" >> $P
export CFM_D_INTAKE=ball CFM_D_COPIES=64 CFM_D_POLAR=through CFM_D_BEAM=0.222 CFM_F_RD=6 CFM_C_SEED=49344
export CFM_D_FIXT=150 CFM_D_CALWARMT=150 CFM_D_CALT=100 CFM_D_LOG=167 CFM_C_DT=0.06 CFM_C_LAMBDA=0.3
export CFM_P_MIX=1,1,1,1 CFM_P_INNER=random CFM_P_SIZE=1 CFM_P_PREFILL=1 CFM_P_CELLS=1
case $RUNG in
  1) RES="CFM_C_H=0.25 CFM_C_NA=8000" ;;
  2) RES="CFM_C_H=0.177 CFM_C_NA=12000" ;;
  3) RES="CFM_C_H=0.125 CFM_C_NA=16000" ;;
esac
PN="CFM_F_PAIR=pn CFM_D_FREEZE=1 CFM_D_STEPS=1667 CFM_D_SEP0=2.05"
PE="CFM_F_PAIR=pe CFM_D_LAYOUT=custom CFM_D_POS=0,0,0;0,0,2.5 CFM_D_MASS=7e6 CFM_D_STEPS=10000"
run() { n=$1; shift; ( export "$@"; $B --ignored --nocapture --exact gpu::dsmc::shell_gates::report_deuteron > $O/$n.txt 2>&1; echo "$n done $(date) exit $?" >> $P ); }
q() { while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do wait -n; done; run "$@" & }
q r${RUNG}_pn_size      $RES $PN
q r${RUNG}_pn_mass_size $RES $PN CFM_P_MASS=1
if [ "$RUNG" = 1 ]; then
  for s in 49344 1 2; do q r1_pe_mass_size_s$s $RES $PE CFM_P_MASS=1 CFM_C_SEED=$s; done
fi
wait
echo "ALL DONE rung $RUNG $(date)" >> $P
