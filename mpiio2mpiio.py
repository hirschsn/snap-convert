#!/usr/bin/python3
#
# mpiio2mpiio.py - Change the number of processes in an MPI-IO snapshot
#
# Copyright 2020 Steffen Hirschmann
#
# Permission to use, copy, modify, and/or distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
# 
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH
# REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY
# AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT,
# INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM
# LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
# OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
# PERFORMANCE OF THIS SOFTWARE.
#

# Description:
# This program converts the offset files of Espresso MPI-IO snapshots to
# a different amount of processes. Subsequently, the newly created snapshot
# is readable by the new amount of processes.
# The program only converts *.pref and *.boff files. All other files are
# re-used.
# This script conveniently creates symbolic links for all files that are
# shared between the new and the old snapshot.
#
# Note: The number of particles is evenly divided between the new number of
# processes. Also, the particles will be out-of-bounds. Use a global resort
# after reading in a snapshot created with this utility.
#
# Usage: mpiio2mpiio OLD_SNAPSHOT_PREFIX NEW_NPROC NEW_SNAPSHOT_PREFIX
# - OLD_SNAPSHOT_PREFIX (string) prefix of the old snapshot (files must exist)
# - NEW_NPROC (int) number of processes the new snapshot should be created for
# - NEW_SNAPSHOT_PREFIX (string) prefix for the newly created files (must not exist yet)
#
# Version: 1.0
#
#

import os
import sys
import numpy as np

def is_sorted(a):
        return np.all(a[1:] >= a[:-1])


if len(sys.argv) != 4 or (len(sys.argv) > 1 and sys.argv[1] == "-h"):
        print("Usage: {} PREFIX NPROC NEW_PREFIX".format(sys.argv[0]), file=sys.stderr)
        raise SystemExit(1)

#PREFIX = "/mnt/neon-scratch/convert-mpiio-hawk/case600_1e-4_norescale-0002080000"
#NEW_PREFIX = "/mnt/neon-scratch/convert-mpiio-hawk/case600_1e-4_norescale-0002080000-256proc"
#TO_PROC = 256

PREFIX = sys.argv[1]
NEW_PREFIX = sys.argv[3]
TO_PROC = int(sys.argv[2])

for suf in ["pref", "boff", "id", "pos", "vel", "type", "head", "bond"]:
        assert(os.path.exists(PREFIX + "." + suf))
        assert(not os.path.exists(NEW_PREFIX + "." + suf))


nglobalpart = os.stat(PREFIX + ".id").st_size // 4 # int32
nproc = os.stat(PREFIX + ".pref").st_size // 4 # int32

if nproc == TO_PROC:
        print("Input snapshot is already for", nproc, "proc", file=sys.stderr)
        raise SystemExit(0)

#
# Load prefixes and calculate the number of particles per process
#
pref = np.fromfile(PREFIX + ".pref", dtype=np.int32)

# Calculate npart per proc (reverse prefix sum)
ppp = np.zeros(pref.shape, dtype=np.int32)
ppp[:-1] = pref[1:] - pref[:-1]
ppp[-1] = nglobalpart - pref[-1]

assert(np.all(ppp > 0))
assert(np.sum(ppp) == nglobalpart)


#
# Bond offsets
#
boff = np.fromfile(PREFIX + ".boff", dtype=np.int32)

assert(boff.shape[0] == nglobalpart + nproc)

# Calculate number of bonds per particle
nbonds_pp = np.zeros(shape=(nglobalpart,), dtype=np.int32)
start = 0
end = 0
for i in range(nproc):
        start = end
        end += ppp[i] + 1
        # Bond offsets are process-local. An easy way to check if we hit the
        # right boundary is to check if the number decreases.
        #print("Process {} start {} end {}".format(i, start, end))
        if i < nproc - 1 and boff[end-1] > 0:
                assert(boff[end-1] > boff[end])
        boff_i = boff[start:end]
        assert(is_sorted(boff_i))
        # Reverse prefix
        bpp_i = np.zeros(boff_i.shape, dtype=np.int32)
        bpp_i = boff_i[1:] - boff_i[:-1]
        
        # Put the number of bonds per particle into the correct spot in "nbonds_pp"
        first_part_id = pref[i]
        last_part_id = pref[i+1] if i < nproc - 1 else nglobalpart
        assert(np.all(nbonds_pp[first_part_id:last_part_id] == 0))
        nbonds_pp[first_part_id:last_part_id] = bpp_i


#
# New prefixes and part per proc
#
new_pref = np.linspace(0, nglobalpart, num=TO_PROC+1, dtype=np.int32)
# This calculation has slightly different indices than the one above because
# we define new_pref to be of size TO_PROC+1 with new_pref[TO_PROC] = nglobalpart.
new_ppp = np.zeros(shape=(new_pref.shape[0] - 1,), dtype=np.int32)
new_ppp = new_pref[1:] - new_pref[:-1]

assert(np.all(new_ppp > 0))
assert(np.sum(new_ppp) == nglobalpart)


#
# Calculate new bond prefixes from "new_ppp" and "nbonds_pp".
#
del boff # Free some RAM
new_boff = np.zeros(shape=(nglobalpart+TO_PROC,), dtype=np.int32)

start = 0
end = 0
for i in range(TO_PROC):
        first_part_id = new_pref[i]
        last_part_id = new_pref[i+1] if i < TO_PROC - 1 else nglobalpart
        assert(last_part_id - first_part_id == new_ppp[i])
        bpp_i = nbonds_pp[first_part_id:last_part_id]
        boff_i = np.zeros(shape=(bpp_i.shape[0]+1,), dtype=np.int32)
        # boff_i[0] := 0
        boff_i[1:] = np.cumsum(bpp_i)
        
        start = end
        end += new_ppp[i] + 1
        assert(end <= new_boff.shape[0])
        new_boff[start:end] = boff_i

new_boff.tofile(NEW_PREFIX + ".boff")
new_pref[:-1].tofile(NEW_PREFIX + ".pref")
for suf in ["head", "id", "pos", "vel", "type", "bond"]: 
        if os.path.exists(PREFIX + "." + suf):
                os.system("ln {}.{} {}.{}".format(PREFIX, suf, NEW_PREFIX, suf))

