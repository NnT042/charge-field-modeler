"""Write one ISABELLE job file: mkjob.py <name> <minutes> <test> KEY=VAL ...
Creates queue/<name>.cmd next to this script (CRLF), for C:\\cfm\\queue."""
import os
import sys

name, minutes, test = sys.argv[1], int(sys.argv[2]), sys.argv[3]
here = os.path.dirname(os.path.abspath(__file__))
os.makedirs(os.path.join(here, 'queue'), exist_ok=True)
lines = ['@echo off', f'rem minutes={minutes}', 'set CFM_PAUSE_FILE=C:\\cfm\\pause\\cfm_pause']
for kv in sys.argv[4:]:
    k, v = kv.split('=', 1)
    assert not any(c in v for c in '&%|<>^"'), v
    lines.append(f'set {k}={v}')
lines.append(f'C:\\cfm\\cfm.exe {test} --ignored --nocapture > C:\\cfm\\out\\{name}.txt 2>&1')
lines.append('exit /b %ERRORLEVEL%')
with open(os.path.join(here, 'queue', name + '.cmd'), 'w', newline='') as f:
    f.write('\r\n'.join(lines) + '\r\n')
print('wrote', name)
