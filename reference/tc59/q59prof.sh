#!/bin/bash
# Doc 15.59 step 0 (docs/VAST_HITBOX_PLAN.md): is the hit-box pile-up the escape boundary?
# Density by radius over time in the escape ball vs the periodic box. h 0.25, n_a 8000, lambda 0.3.
export PATH=/usr/bin:/bin:$PATH
cd /c/Users/nnt04/Project/charge-field-modeler
P=reference/tc59/progress_prof.txt
run() { name=$1; shift; env "$@" cargo test --manifest-path rust/Cargo.toml --profile fast --features gpu report_photon_mix_profile -- --ignored --nocapture --test-threads=1 2>&1 | grep -E "===|^t |wall|panicked|periodic fill|columns|MAX_CELL" > reference/tc59/prof_$name.txt; echo "$name done $(date)" >> $P; }
echo "START $(date)" > $P
EQ="CFM_P_MIX=1,1,1,1 CFM_P_INNER=random"
run esc_k0 CFM_P_TIME=120
run esc_size $EQ CFM_P_SIZE=1 CFM_P_TIME=120
run per_size $EQ CFM_P_SIZE=1 CFM_P_BOX=periodic CFM_P_TIME=60
run per_mass_size $EQ CFM_P_SIZE=1 CFM_P_MASS=1 CFM_P_BOX=periodic CFM_P_TIME=30 CFM_P_SNAP=5
run esc_mass_size $EQ CFM_P_SIZE=1 CFM_P_MASS=1 CFM_P_TIME=60 CFM_P_SNAP=5
echo "ALL DONE $(date)" >> $P
