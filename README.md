# Binary Data Layout and Interpretation

Steffen Hirschmann

01.03.2021

The simulation data of the 109.85 mio. particle simulation is 
stored in a custom binary, parallel format as created by 
ESPResSo's MPI-IO routine. The documentation of the ESPResSo 
functionality to generate and read it can be found here: 
http://espressomd.org/html/doc/io.html#writing-mpi-io-binary-files

The purpose of this document is to clarify the interpretation and 
usage of this kind of data.

## Data format

Each MPI-IO snapshot consists of at least 3 and at most 8 files. 
The files correspond to different data fields from the particles. 
A full snapshot with of particle positions, velocities, types, 
and bonding information consists of the following files:

- snapshot-1.head – Administrative information (stores the data 
  fields that the snapshot contains).
- snapshot-1.pref – Administrative information to read the 
  snapshot in parallel.
- snapshot-1.id – IDs of all particles. (Simple array, int32_t)
- snapshot-1.pos – Positions of all particles. (Simple array, 
  double[3])
- snapshot-1.vel – Velocities of all particles. (Simple array, 
  int32_t)
- snapshot-1.type – The particle types. (Simple array, double[3])
- snapshot-1.boff – Administrative information: Number of bonds 
  per particle (Not a simple array)
- snapshot-1.bond – Bonding information (Not a simple array)

The files termed “simple array” above can be interpreted like:

```
int type = type_file[id_file[myid]]; // (1)
double *pos = &pos_field[3 * id_file[myid]]; // (2)
printf("Particle %i, type: %i, pos: %lf %lf %lf\n", myid, type, pos[0], pos[1], pos[2]);
```

In this example, line (1) reads the type of particle “myid” given 
that “type_field” is the data from file “snapshot-1.type” read as 
an array of integer, see e.g. C-function fread (linux shell: `man 
3 fread`). Line (2) shows how to get position/velocity 
information, which is 3 doubles per particle.

The data labelled “not a simple array” can not be read in like 
shown above because it depends on the contents of “
snapshot-1.pref”, which, in turn, stores a prefix sum of the 
number of particles per process. Parsing these files is better 
left to the tools described below. Feel free to copy the parsing 
code from any of these tools. Possibly the best option is to copy 
snapshot.hpp from snap-analyze (
https://github.com/hirschsn/snap-analyze/blob/master/util/snapshot.hpp
, see the #include directives on which other header files are 
needed) and to use the template function “snapshot_iter”.

## Tools

### Mpiio2blockfile.py

Mpiio2blockfile.py (https://github.com/hirschsn/snap-convert) is 
a Python2(!) script that converts MPI-IO into textual data that 
ESPResSo used before version 4.0. It cannot be read by ESPResSo 
anymore. The data can, however, be read by humans. The data that 
this script outputs looks like this:

```
{particles {id type pos v}
	{0 0 10.0 20.0 30.0 -1.0 -2.0 2.1}
	{1 0 12.3 0.5 7.3 -0.26 -0.5 -0.6}
	...
}
{bonds
	{0 { {0 42965} {151 35837 42965} {0 42560} {65 42965 42560} } }
	{1 { {0 33682} {0 521} {121 41644 521} } }
}
```


The first part is the particle information and the second the 
bond information. Each line corresponds to one particle and each 
line contains data for “id type pos v”, in this order. Thus, the 
first line is the particle with id 0, its type is 0, its position 
is (10, 20, 30) and its velocity is (-1, -2, 2.1), where vectors 
are ordered (x, y, z).

Each line of the bonding information stores the bonds of one 
particle. The first line contains a list of bonds of the 0-th 
particle: Bond type 0 to particle with id 42965, bond type 151 (a 
three particle bond) with particle ids 35837 and 42965 and so on. 
It is important to note that ESPResSo places bonds between 
particles a and b only on one of the two particles. Thus, if you 
find (bond-type, b) in a's bond list, there won't be a 
corresponding (bond-type, a) entry in b's bond list.



### mpiio2mpiio.py

If an MPI-IO snapshot was dumped with an ESPResSo running on 128 
processes, it can only be read in, again, by ESPResSo also 
running on 128 processes. Mpiio2mpiio.py (https://github.com/hirschsn/snap-convert
) is a Python3(!) script that can change this. Called like

```
mpiio2mpiio.py snapshot-1 8192 snapshot-1-8192proc
```

It generates a new MPI-IO snapshot called “snapshot-1-8192proc” 
from “snapshot-1” that can be read in by an ESPResSo instance 
running on 8192 processes (second argument to the script 
specifies the desired number of processes).

The script changes the order of the particles and almost all of 
them will be out of box on the process that reads them. This 
causes no problem because ESPResSo performs a global particle 
exchange on the first timestep, should this occur.

### snap-analyze

Snap-analyze (https://github.com/hirschsn/snap-analyze) is a 
program written in C++ that analyzes the bonding structure of an 
MPI-IO snapshot. It has two modes: The first (`--df`, requires also: `--sigma 1.0`
 and `--box 2600.0`; insert the proper values for the box length and sigma) 
extracts all aggregates from a snapshot, calculates their radius 
of gyration and fractal dimension. The result is printed to 
stdout.

The second mode generates a file for each aggregate that contains 
the particle positions of the particles that are part of the 
aggregate. Each individual file can be read in, e.g. with numpy's 
loadtxt. Note, that this potentially generates a HUGE number of 
files. For a snapshot from the 109 mio particle simulation, the 
number of files will be too long for the Linux command line: E.g. 
`rm *` won't work because the argument list gets too long. Here,
`find . -maxdepth 1 -type f -iname "*_POS_*" -exec rm {} \;`
is your fried. Substitute whatever command (instead “rm”) you 
want to run. “{}” is replaced by “find” with each individual file 
name.
