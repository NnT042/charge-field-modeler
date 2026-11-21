#!/bin/bash
# Pause (stop) this Vast instance when the queue is finished.
# Trigger: "ALL DONE" in ~/out/progress.txt, or the queue's tmux session 'q' has gone away,
# or a hard cap of MAX_H hours. Then wait GRACE_MIN minutes (so ESTELLE can copy results)
# and stop. Stop keeps the disk; destroy stays a manual step.
# Start:  tmux new -d -s stop 'bash ~/autostop.sh'
GRACE_MIN=${GRACE_MIN:-30}; MAX_H=${MAX_H:-12}
L=$HOME/out/autostop.log; mkdir -p $HOME/out
t0=$(date +%s)
echo "autostop armed $(date), grace ${GRACE_MIN} min, cap ${MAX_H} h" >> $L
while true; do
  why=""
  grep -q "ALL DONE" $HOME/out/progress.txt 2>/dev/null && why="queue done"
  [ -z "$why" ] && ! tmux has-session -t q 2>/dev/null && why="queue session gone"
  [ -z "$why" ] && [ $(( $(date +%s) - t0 )) -gt $(( MAX_H * 3600 )) ] && why="hard cap"
  [ -n "$why" ] && break
  sleep 60
done
echo "trigger: $why at $(date); stopping in ${GRACE_MIN} min" >> $L
sleep $(( GRACE_MIN * 60 ))
echo "stopping $(date)" >> $L
vastai stop instance "$CONTAINER_ID" --api-key "$CONTAINER_API_KEY" >> $L 2>&1
