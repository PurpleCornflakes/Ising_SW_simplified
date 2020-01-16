# Potts Model with Swendsen-Wang Algorithm

## Potts Model

## Swendsen-Wang Cluster Algorithm

## Much More efficient than Metropolis Algorithm

## Code
1. iterate SW to a equilibrium state, meanwhile record *energy* and *magnetization* at each step, check their histogram, output a spin lattice in equilibrium to a file.
2. read in the spin-lattice file, iterate and record slices at each step. 
3. draw MutualInfo versus distance.

gcc -o sw -lm SW_Potts2D.c

qsub T*0.pbs
# ./sw infile0* 0
qsub T*1.pbs
# ./sw infile1* 1

infile* include arguments for main().

e_dist* is energy recorded.
m_dist* is magnetic moment recorded.

slices* is the middle row of the spin_lattice.
spin_lattice* is the last state of the spin lattice.
