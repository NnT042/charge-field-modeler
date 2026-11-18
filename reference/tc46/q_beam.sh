#!/bin/bash
# Doc 15.46 overnight queue: neutron trapped-release beam width vs the trace's clear polar hole.
# Measured (report_trace_clear_hole): ball centres never come closer than 0.347 to the axis, the
# ball is 0.125, so the hole the trace leaves clear is 0.222. The old 1/3 beam = the centres' hole.
# Trace bodies: 64 copies, surface emission, through release, lambda 0.3 (same as tc45 q_through).
export PATH=/usr/bin:/bin:$PATH
cd /c/Users/nnt04/Project/charge-field-modeler/rust
B=target/fast/deps/charge_field_modeler-8ba177252dd33c00.exe
O=../reference/tc46
P=$O/progress.txt
export CFM_D_INTAKE=ball CFM_D_COPIES=64 CFM_D_POLAR=through CFM_C_LAMBDA=0.3 CFM_C_H=0.25 CFM_C_DT=0.06 CFM_F_RD=6
frozen() { export CFM_D_FREEZE=1 CFM_D_FIXT=150 CFM_D_CALWARMT=150 CFM_D_CALT=100 CFM_D_STEPS=1667 CFM_D_LOG=167; }
run() { $B --ignored --nocapture --exact gpu::dsmc::shell_gates::report_deuteron > $O/$1.txt 2>&1; echo "$1 done $(date)" >> $P; }
fself() {
  F=$(grep "^self-force neutron" $O/$1.txt | head -1 | sed -E 's/.*\(([^)]*)\).*/\1/; s/[ +]//g')
  T=$(grep "^self-torque neutron" $O/$1.txt | head -1 | sed -E 's/.*\(([^)]*)\).*/\1/; s/[ +]//g')
  echo "proton:-106.9,-12.0,22.6,1.5,1.8,-4.2;neutron:$F,$T"
}
echo "START $(date)" > $P
frozen
# 1. Deuteron at 3 with the 0.222 beam (also calibrates the neutron self-force at this beam).
CFM_D_BEAM=0.222 CFM_F_PAIR=pn CFM_D_SEP0=3 run pn_axis_3_b222
FS222=$(fself pn_axis_3_b222); echo "FSELF b222 $FS222" >> $P
# 2. He3 at contact, 0.222 beam, three seeds; matched 1/3-beam controls on the two new seeds
#    (seed-0 1/3 result is tc45 tr_he3_axis_205_through: -1615/-1553).
FS333=$(grep -h "^FSELF" ../reference/tc45/progress_through.txt | head -1 | sed 's/^FSELF //')
echo "FSELF b333 (tc45) $FS333" >> $P
CFM_D_BEAM=0.222 CFM_F_PAIR=pnp CFM_D_SEPS=2.05,2.05 CFM_D_FSELF="$FS222" run he3_205_b222_s0
CFM_C_SEED=777   CFM_D_BEAM=0.222 CFM_F_PAIR=pnp CFM_D_SEPS=2.05,2.05 CFM_D_FSELF="$FS222" run he3_205_b222_s777
CFM_C_SEED=777   CFM_D_BEAM=0.3333333 CFM_F_PAIR=pnp CFM_D_SEPS=2.05,2.05 CFM_D_FSELF="$FS333" run he3_205_b333_s777
CFM_C_SEED=31337 CFM_D_BEAM=0.222 CFM_F_PAIR=pnp CFM_D_SEPS=2.05,2.05 CFM_D_FSELF="$FS222" run he3_205_b222_s31337
CFM_C_SEED=31337 CFM_D_BEAM=0.3333333 CFM_F_PAIR=pnp CFM_D_SEPS=2.05,2.05 CFM_D_FSELF="$FS333" run he3_205_b333_s31337
# 3. Dose response: a narrower beam (0.15) on seed 0.
CFM_D_BEAM=0.15 CFM_F_PAIR=pnp CFM_D_SEPS=2.05,2.05 CFM_D_FSELF="$FS222" run he3_205_b150_s0
# 4. Free He3 from contact with the 0.222 beam: 250 tu held, 600 tu free (ch46 masses).
unset CFM_D_FREEZE CFM_D_CALWARMT CFM_D_CALT
export CFM_D_FIXT=250 CFM_D_STEPS=14167 CFM_D_LOG=500 CFM_D_LAYOUT=custom
export CFM_D_MASS=7e6 CFM_D_INERTIA=2.8e6 CFM_D_LSPIN=8.75e5
CFM_D_BEAM=0.222 CFM_F_PAIR=pnp CFM_D_POS="0,0,-2.05;0,0,0;0,0,2.05" CFM_D_FSELF="$FS222" run he3_free_b222_s0
echo "ALL DONE $(date)" >> $P
