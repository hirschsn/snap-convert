#!/usr/bin/python2
#
# mpiio2blockfile.py - Transform an MPI-IO snapshot dumped by ESPResSo into a
# human readable format.
#
#
# Copyright (C) 2016-2017 Steffen Hirschmann
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
#
from __future__ import print_function

import os
import sys
import array

if len(sys.argv) < 2:
    sys.stderr.write("Usage: {} PREFIX\n".format(sys.argv[0]))
    sys.exit(1)

fpref = sys.argv[1]
headf = fpref + ".head"
preff = fpref + ".pref"
idf = fpref + ".id"
typef = fpref + ".type"
posf = fpref + ".pos"
velf = fpref + ".vel"
bofff = fpref + ".boff"
bondf = fpref + ".bond"

read_typ = os.access(typef, os.R_OK)
read_pos = os.access(posf, os.R_OK)
read_vel = os.access(velf, os.R_OK)
read_bonds = os.access(bofff, os.R_OK)

# Determine nproc, etc. at time of writing
nproc = os.stat(preff).st_size / 4
ntotalpart = os.stat(idf).st_size / 4
if read_bonds:
    ntotalbond = os.stat(bondf).st_size / 4

# Read header - fields, n_bonded_ia, bonded_ia_params[:].nums
fields = 0
nbia = 0
biaparams = array.array("i")
with open(headf) as f:
    afields = array.array("I")
    afields.read(f, 1)
    fields = afields[0]
    anbia = array.array("i")
    anbia.read(f, 1)
    nbia = anbia[0]
    biaparams.read(f, nbia)

# Read prefixes
pref = array.array("i")
with open(preff) as f:
    pref.read(f, nproc)

# Read ids
idn = array.array("i")
with open(idf) as f:
    idn.read(f, ntotalpart)

# Read types
if read_typ:
    typ = array.array("i")
    with open(typef) as f:
        typ.read(f, ntotalpart)

# Read pos
if read_pos:
    pos = array.array("d")
    with open(posf) as f:
        pos.read(f, 3 * ntotalpart)

# Read vel
if read_vel:
    vel = array.array("d")
    with open(velf) as f:
        vel.read(f, 3 * ntotalpart)

# Read bonds
if read_bonds:
    boff = array.array("i")
    with open(bofff) as f:
        boff.read(f, ntotalpart + nproc)

    bond = array.array("i")
    with open(bondf) as f:
        bond.read(f, ntotalbond)


# Print particles in blockfile format
spec = ["id"]
if read_typ:
    spec.append("type")
if read_pos:
    spec.append("pos")
if read_vel:
    spec.append("v")

print("{particles {%s}" % " ".join(spec))
for i in range(ntotalpart):
    print("\t{%i" % idn[i], end="")
    if read_typ:
        print(" %i" % typ[i], end="")
    if read_pos:
        print(" %r %r %r" % (pos[3*i], pos[3*i+1], pos[3*i+2]), end="")
    if read_vel:
        print(" %r %r %r" % (vel[3*i], vel[3*i+1], vel[3*i+2]), end="")
    print("}")
print("}")

if not read_bonds:
    sys.exit(0)

# Print bonds in blockfile format
print("{bonds")
addend = 0 # ntotal bonds of previous processors
for rank in range(nproc):
    # The start and end indices for the boff array are determined via
    # pref. However, there are (nlocalpart + 1) boff entries per proc.
    start = pref[rank] + rank
    end = rank + (pref[rank + 1] if rank < nproc - 1 else ntotalpart)
    for pid, i in enumerate(range(start, end)):
        print("\t{%i { " % idn[pref[rank] + pid], end="")
        # The start and end indices for the bond array are determined
        # via boff. However, boff does only *locally* store prefixes,
        # i.e. they have to be globalized by adding the total number of
        # bonds on all ranks before this one.
        j = addend + boff[i]
        while j < addend + boff[i + 1]:
            bond_num = bond[j]
            j += 1
            npartners = biaparams[bond_num]
            print("{%i" % bond_num, end="")
            for _ in range(npartners):
                print(" %i" % bond[j], end="")
                j += 1
            print("} ", end="")
        print("} }")
    addend += boff[end]
print("}")
