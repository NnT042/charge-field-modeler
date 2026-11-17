#!/bin/bash
# Doc 15.45 gate 1: lone proton and lone neutron as trace bodies (Intake::Ball, maze, surface
# emission), opacity = trace copies 1/4/16/64. lambda 0.3, n_a 8000, h 0.25, dt 0.06, R_D 6,
# 100 tu warm-up + 100 tu averaged.
cd /c/Users/nnt04/Project/charge-field-modeler/rust
B=target/fast/deps/charge_field_modeler-8ba177252dd33c00.exe
O=../reference/tc45
export CFM_I_RUNS="ball:swept:0" CFM_C_LAMBDA=0.3 CFM_C_NA=8000 CFM_C_H=0.25 CFM_C_DT=0.06 CFM_C_WARMUP=1667 CFM_C_AVG=1667 CFM_I_SPEED=surface
for cl in proton neutron; do for c in 1 4 16 64; do
  CFM_I_CLASS=$cl CFM_C_COPIES=$c $B --ignored --nocapture --exact gpu::dsmc::shell_gates::report_intake > $O/lone_${cl}_c$c.txt 2>&1
  echo "$cl c$c done $(date)" >> $O/progress.txt
done; done
echo "ALL DONE $(date)" >> $O/progress.txt
