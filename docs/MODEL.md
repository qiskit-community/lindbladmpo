# Lindblad dynamics

The state of an open quantum system is defined by a density matrix <img src="https://render.githubusercontent.com/render/math?math=\rho" style="vertical-align:bottom">.
In the case where the system is made of qubits (two-level systems), the Hilbert space has dimension 2^N (and the density matrix is 2^N by 2^N in size).

The time-evolution of the system coupled to a Markovian environment (a memory-less bath) is often described using the Lindblad master equation,

<img src="https://render.githubusercontent.com/render/math?math=\frac{\partial}{\partial t}\rho = -\frac{i}{\hbar}[H,\rho]%2b\mathcal{D}[\rho]">

where <img src="https://render.githubusercontent.com/render/math?math=[\cdot,\cdot]" style="vertical-align:bottom"> denotes the commutator of two operators.

The first term in the equation above generates the unitary evolution due to interactions and coherent terms in the Hamiltonian, while the second term (sometimes known as the dissipator, or Lindbladian) is a superoperator that accounts for incoherent dephasing and relaxation processes due to the environment.

# The simulator

The solver integrates numerically the dynamics in fixed time steps <img src="https://render.githubusercontent.com/render/math?math=\tau">, to obtain a solution to the master equation as a function of time. A Trotter expansion of (up-to) the fourth order in the time step is being employed.
The solver is based on matrix product states (MPS) and matrix product operators (MPO).

## The Model

The solver can model interacting qubits in a lattice with an arbitrary connectivity defined by the user.

The simulator uses a fixed Hamiltonian with coefficients which can be defined by the user. The Hamiltonian is defined in the rotating frame with
respect to the frequency of an drive of a uniform frequency applied to all qubits (which are not necessarily identical themselves, however, and may be driven at different amplitudes and phases), or not driven at all.
The Hamiltonian is represented as the sum of on-site terms and an interaction part,

<img src="https://render.githubusercontent.com/render/math?math={H}/{\hbar} = \sum_{i}\frac{1}{2}\left[h_{z,i}\sigma_i^z  %2b h_{x,i}\sigma_i^x %2b h_{y,i}\sigma_i^y\right] %2b T,">

with the interaction being,

<img src="https://render.githubusercontent.com/render/math?math=T = \sum_{ i}^N\sum_{ j\neq i}^N \left(J_{ij}\sigma^%2b_i \sigma^-_{j} %2b {\rm h.c.} %2b \frac{1}{2} J_{ij}^z \sigma^z_i \sigma^z_{j}\right)=\\ \frac{1}{2}\sum_{ i}^N\sum_{ j\neq i}^N \left(J_{ij}\sigma^x_i \sigma^x_{j} %2b J_{ij}\sigma^y_i \sigma^y_{j} %2b J_{ij}^z \sigma^z_i \sigma^z_{j}\right),">

where

<img src="https://render.githubusercontent.com/render/math?math=\sigma_i^a, a=\{x,y,z\}," style="vertical-align:bottom">

are the Pauli matrices at each site and

<img src="https://render.githubusercontent.com/render/math?math=\sigma^{\pm}_i = {\sigma^{x}_i\pm i\sigma^y_i}/{2}," style="vertical-align:bottom">

are the ladder operators.

The Lindbladian used in the simulator is also fixed with coefficients that can be defined by the user. It is represented as a sum of three terms,

<img src="https://render.githubusercontent.com/render/math?math=\mathcal{D} = \sum_j \mathcal{D}_j," style="vertical-align:bottom">

The first two terms model exchange with a thermal bath:

<img src="https://render.githubusercontent.com/render/math?math=\mathcal{D}_0[\rho] = \sum_i g_{0,i}\left(\sigma_i^%2b \rho\sigma_i^- - \frac{1}{2} \{\sigma_i^- \sigma_i^%2b,\rho\}\right),">

<img src="https://render.githubusercontent.com/render/math?math=\mathcal{D}_1[\rho]=\sum_i g_{1,i}\left( \sigma_i^-\rho \sigma_i^{%2b}-\frac{1}{2}\left\{\sigma_i^{%2b}\sigma_i^-,\rho\right\}\right),">

and the last term models pure dephasing in xy plane,

<img src="https://render.githubusercontent.com/render/math?math=\mathcal{D}_2[\rho] = \sum_i g_{2,i} \left(\sigma_i^z \rho\sigma_i^z - \rho\right).">
