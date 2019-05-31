Name of Model: Nodes of Ranvier with Left shifted Nav channels

Reference: Boucher PA, Joos B, Morris CE (2012) Coupled left-shift of
Nav channels: modeling the Na(+)-loading and dysfunctional
excitability of damaged axons J Comput Neurosci.

ModelDB Accession: 144481
These files were supplied by Dr. B. Joos.

The two programs CLSRanvier.f and propagation.f simulate the
excitability of a myelinated axon with injured nodes of Ranvier.  The
injury is simulated as the Coupled Left shift (CLS) of the activation
and inactivation of a fraction of Nav channels.

Program CLSRanvier.f

This program calculates the time evolution of the membrane potential
at a node of Ranvier with left-shifted Nav currents.

Program propagation.f

This program calculates the time evolution of the membrane potential
for a series of nodes of Ranvier connected by simple resistors,
representing the myelenated segments of axons.  The simulation is
carried 30 times (30 values of shift) for 6 different values of LS. No
stimulation current is used on the initial node, so spontaneous firing
of the nodes is studied.  A single node (node 5) has shifted channels.

Usage:

These programs are compatible with the linux gfortran compiler.  They
can be compiled with the commands

gfortran -o CLSRanvier CLSRanvier.f
gfortran -o propagation propagation.f

and then run with commands

./CLSRanvier
./propagation

On a 2.40GHz Intel Core 2 Duo P8600 (Dell laptop) it takes 15 seconds
to run CLSRanvier and 11 minutes to run propagation.
