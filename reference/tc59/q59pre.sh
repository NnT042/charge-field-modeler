#!/bin/bash
# Doc 15.59: pre-filled start (CFM_P_PREFILL=1). Profile check of the worst case, then the size-only p-n pair.
# ISABELLE runs q59_d (size + mass p-n, pre-filled).
export PATH=/usr/bin:/bin:$PATH
cd /c/Users/nnt04/Project/charge-field-modeler
P=reference/tc59/progress_pre.txt
echo "START $(date)" > $P
EQ="CFM_P_MIX=1,1,1,1 CFM_P_INNER=random"
env $EQ CFM_P_SIZE=1 CFM_P_MASS=1 CFM_P_PREFILL=1 CFM_P_TIME=30 CFM_P_SNAP=5 cargo test --manifest-path rust/Cargo.toml --profile fast --features gpu report_photon_mix_profile -- --ignored --nocapture --test-threads=1 2>&1 | grep -E "===|^t |wall|panicked|prefill|MAX_CELL" > reference/tc59/prof_esc_mass_size_prefill.txt
echo "prof done $(date)" >> $P
cd rust
B=/c/Users/nnt04/AppData/Local/Temp/claude/63a7c279081ac85c/b395b8e1-54e5-4c45-a8ec-a703a023d0dd/scratchpad/cfm59b.exe
export CFM_D_INTAKE=ball CFM_D_COPIES=64 CFM_D_POLAR=through CFM_D_BEAM=0.222 CFM_F_RD=6 CFM_C_SEED=49344
export CFM_D_FIXT=150 CFM_D_CALWARMT=150 CFM_D_CALT=100 CFM_D_LOG=167 CFM_C_DT=0.06 CFM_C_LAMBDA=0.3 CFM_C_H=0.25
env $EQ CFM_P_SIZE=1 CFM_P_PREFILL=1 CFM_F_PAIR=pn CFM_D_FREEZE=1 CFM_D_STEPS=1667 CFM_D_SEP0=2.05 $B --ignored --nocapture --exact gpu::dsmc::shell_gates::report_deuteron > ../reference/tc59/pn_eq_random_size_prefill.txt 2>&1
echo "pn size prefill done $(date)" >> ../$P
echo "ALL DONE $(date)" >> ../$P
