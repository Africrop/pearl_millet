//Number of population samples (demes)
4 number of populations
// Demes sizes
NCult$
Nw$
Ne$
Nc$
//Sample sizes (followed by sampling times if required)
9
5
4
9
//Growth rates : negative growths implies population expansion
GRate$
0
0
0
//Number of migration matrices: 0 implies no migration between demes
0
//Historical events: first line number of hist events, then each line is an event with time, source, sink, migrants, new size, growth rates, migr.matrix
3 historical event
T0$ 0 3 1 1 0 0
T1$ 1 3 1 Res1a$ 0 0
T2$ 3 2 1 Res2a$ 0 0
//Number of independant loci (chromosomes)
1 0
//Per chromosome: Number of linkage blocks
1
//per block: data type, num loci, rec.rate and mut rate + optional parameters
FREQ 1 0 1.7e-8
