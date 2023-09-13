# 2D Ising model
Parallel implementation of the 2D Ising model in C

---
## Background

The (nearest neighbour) Ising model, in absence of external magnetic fields, is described by the following Hamiltonian
$$\mathcal{H}=-J\sum _{&lt; i,j &gt;} s_is_j.$$

Here $s_i\in$ \{ $+1,-1$ \} is the $i$-th spin on the underlying, usually square, lattice, and &lt; $i,j$ &gt; denotes a sum over nearest neighbours. To simplify the model, only nearest neighbour interactions are considered, the strength of which is described by $J$, which is positive for ferromagnetic coupling. It is one of the most studied and well-known models in statistical physics. The 1D Ising model, which is a chain of connected spins, does not allow for a transition from a ferromagnetic- to a paramagnetic phase, since the ferromagnetic state is not the ground state, i.e. the state with the lowest energy at $T=0$. However, the ferromagnetic state is the ground state for the 2D Ising model and is subject to such a phase transition. At non-zero temperatures, spins randomly flip, which naturally happens more often for increasing temperatures, until a certain temperature $T_C$ is reached, after which the order completely vanishes. Typical for such statistical systems is the presence of an order parameter, which for the Ising model is given by the total magnetization of the system, $M$, and is hence of high importance. This critical temperature can be derived analytically in 2D by explicit computation of the so-called partition function and gives a value of $$T_C=\frac{2J}{\text{ln}(1+\sqrt{2})},$$
 which gives $T_C \approx$ 2.269 for $J=1$.


## Implentation

Here, the dynamics of the Ising model are simulated. The critical temperature is determined by studying the total magnetization at a large range of temperatures and looking at where anomalous behaviour occurs. Simulations of the Ising model are done using the Metropolis algorithm, which mimics the stochastic nature of the spins at non-zero temperatures. Parallelization is done using OpenMP, which divides the lattice into vertical strips, and runs the Metropolis algorithm over each of these strips, individually. 

For a large range of temperatures, the ensemble-averaged magnetization takes the following shape.
<p align="center">
<img src="https://user-images.githubusercontent.com/98324298/226294129-7123b163-3fd8-4940-a7f7-91bd499ae40c.png" width="600">
</p>

For low $T$ a ferromagnetic configuration is indeed observed, which slowly decays until $T_C$ is reached, after which the state is paramagnetic. Magnetizations at temperatures close to $T_C$ fluctuate a bit as there is a small bias when setting up the spins of the system to help reach thermal equilibrium in the initialization Metropolis steps of the simulations. Random numbers are generated using the fast [xoroshiro128+](https://prng.di.unimi.it/) generator [1].


## 
The code is built with the `makefile` that is supplied. The binary works like 

`./ising n_thrm n_measure N_threads`

Here `n_thrm` defines the number of sweeps in the thermalization stage, `n_measure` defines the number of sweeps in the measurement stage, and `N_threads` is the number of cores when ran in parallel. One sweep is defined as having attempted to flip every spin on the lattice once.

## Future improvements
- Divide up the lattice into blocks instead of strips for better cache performance.
- Implement Wolff algorithm to collectively flip spins for more accurate results.
- Parallelize using MPI.

## References
[1] D. Blackman and S. Vigna, “Scrambled linear pseudorandom number generators,”
ACM Transactions on Mathematical Software (TOMS), vol. 47, no. 4, pp. 1–32,
2021.
