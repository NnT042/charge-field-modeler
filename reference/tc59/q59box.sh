#!/bin/bash
# tc59 field-box runs (doc 15.59): no bodies, photon mix knobs.
export PATH=/usr/bin:/bin:$PATH
cd /c/Users/nnt04/Project/charge-field-modeler
run() { name=$1; shift; env "$@" cargo test --manifest-path rust/Cargo.toml --profile fast --features gpu report_photon_mix -- --ignored --nocapture --test-threads=1 2>&1 | grep -E "===|census|per unit|by radius|^  r |wall|panicked" > reference/tc59/box_$name.txt; echo "$name done $(date)" >> reference/tc59/progress.txt; }
echo "START $(date)" > reference/tc59/progress.txt
for inner in match plus minus random; do run eq_${inner} CFM_P_MIX=1,1,1,1 CFM_P_INNER=$inner; done
for inner in plus minus random; do run eq_${inner}_mass_size CFM_P_MIX=1,1,1,1 CFM_P_INNER=$inner CFM_P_MASS=1 CFM_P_SIZE=1; done
run eq_random_mass CFM_P_MIX=1,1,1,1 CFM_P_INNER=random CFM_P_MASS=1
run eq_random_size CFM_P_MIX=1,1,1,1 CFM_P_INNER=random CFM_P_SIZE=1
run k0_size CFM_P_MIX=1,0,0,0 CFM_P_SIZE=1 CFM_P_CEIL=3
echo "ALL DONE $(date)" >> reference/tc59/progress.txt
