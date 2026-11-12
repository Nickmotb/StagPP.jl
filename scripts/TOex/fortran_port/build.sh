#!/bin/bash
# Builds test_toex against a SCRATCH watercarb build in $WC (never the working copy): make stagyympi there first.
#  Run: ./test_toex parity|closure|cell|speed <tables dir>; parity needs states.txt and j_*.txt from ref_coh.jl.
set -e
WC=${WC:-$TMPDIR/wc}
LINK=$(grep -h "\-o stagyympi" $WC/build*.log | head -1)
OBJS=$(echo "$LINK" | tr ' ' '\n' | grep '\.o$' | grep -v '^main.o$' | sed "s|^|$WC/|")
mpifort -O3 -g -cpp -ffree-line-length-none -fdefault-double-8 -fdefault-real-8 -I$WC -c test_toex.f90 -o test_toex.o
mpifort -O3 -g -fdefault-double-8 -fdefault-real-8 -o test_toex test_toex.o $OBJS ${LINK#*gpukernels_mod.o} 2>&1 | grep -v duplicate || true
