# Lindblad dynamics

The state of an open quantum system is defined by a density matrix <img src="https://render.githubusercontent.com/render/math?math=\rho">.
In the case where the system is made of of qubits (two level systems), the Hilbert space has dimension 2^N (therefore the density matrix is 2^N by 2^N in size).

When the system is coupled to the environment (memory-less bath), it's state evolves over time.
The notation of such evolution is <img src="https://render.githubusercontent.com/render/math?math=\mathcal{L}[\rho]" style="vertical-align:middle">.
This Liouvillian is a (linear) superoperator that maps the operator <img src="https://render.githubusercontent.com/render/math?math=\rho" style="vertical-align:middle"> to another operator - the time derivative of <img src="https://render.githubusercontent.com/render/math?math=\rho" style="vertical-align:middle">, and is often described using the master equation:

<img src="https://render.githubusercontent.com/render/math?math=\frac{\partial}{\partial t}\rho = \mathcal{L}[\rho]\equiv -\frac{i}{\hbar}[H,\rho]%2b\mathcal{D}[\rho]">

where <img src="https://render.githubusercontent.com/render/math?math=[\cdot,\cdot]" style="vertical-align:middle"> describes the commutator of two operators.

The first part of the equation, <img src="https://render.githubusercontent.com/render/math?math=\frac{i}{\hbar}[H,\rho]" style="vertical-align:middle">, expresses the unitary evolution due to interactions and coherent driving terms,
while <img src="https://render.githubusercontent.com/render/math?math=\mathcal{D}[\rho]" style="vertical-align:middle">  (sometimes known as the dissipator, or Lindbladian) is a superoperator that accounts for incoherent dephasing and relaxation processes due to the environment.

# The simulator

The simulator simulates numerically Lindbladian dynamics and offers approximate solution to the master equation.
The simulator's solver is based on matrix product states (MPS) and matrix product operators (MPO).
The matrix product state is used to represent the state of the system, while the matrix product operator is used to represent the time-evolution operator.
The solver integrates the dynamics in fixed time steps <img src="https://render.githubusercontent.com/render/math?math=\tau"> defined by the user.

## The Model

The simulator can model interacting qubits in a planar lattice with an arbitrary connectivity defined by the user.

The simulator uses a fixed Hamiltonian with coefficients which can be defined by the user. The Hamiltonian is defined in the rotating frame with
respect to the frequency of an identical drive applied to all qubits (which are not necessarily identical themselves, however, and may be driven at different amplitudes).
The Hamiltonian is represented as the sum of on-site terms and an interaction part T:

<img src="https://render.githubusercontent.com/render/math?math=\frac{H}{\hbar} = \sum_{i}\frac{1}{2}\left[h_{z,i}\sigma_i^z  %2b h_{x,i}\sigma_i^x %2b h_{y,i}\sigma_i^y\right] %2b T">

with the interaction being:

<img src="https://render.githubusercontent.com/render/math?math=T = \sum_{ i}^N\sum_{ j\neq i}^N \left(J_{ij}\sigma^%2b_i \sigma^-_{j} %2b {\rm h.c.} %2b \frac{1}{2} J_{ij}^z \sigma^z_i \sigma^z_{j}\right)=\\ \frac{1}{2}\sum_{ i}^N\sum_{ j\neq i}^N \left(J_{ij}\sigma^x_i \sigma^x_{j} %2b J_{ij}\sigma^y_i \sigma^y_{j} %2b J_{ij}^z \sigma^z_i \sigma^z_{j}\right)">

while <img src="https://render.githubusercontent.com/render/math?math=\sigma_i^a, a=\{x,y,z\}" style="vertical-align:middle"> are the Pauli matrices at each site and <img src="https://render.githubusercontent.com/render/math?math=\sigma^{\pm}_i = {\sigma^{x}_i\pm i\sigma^y_i}/{2}" style="vertical-align:middle"> are the ladder operators.

The Lindbladian used in the simulator is also fixed with coefficients which can be defined by the user. It is represented as a sum of three terms -
<img src="https://render.githubusercontent.com/render/math?math=\mathcal{D} = \sum_j \mathcal{D}_j" style="vertical-align:middle">

The first two terms represent exchange with a thermal bath:

<img src="https://render.githubusercontent.com/render/math?math=\mathcal{D}_0[\rho] = \sum_i g_{0,i}\left(\sigma_i^%2b \rho\sigma_i^- - \frac{1}{2} \{\sigma_i^- \sigma_i^%2b,\rho\}\right)">

<img src="https://render.githubusercontent.com/render/math?math=\mathcal{D}_1[\rho]=\sum_i g_{1,i}\left( \sigma_i^-\rho \sigma_i^{%2b}-\frac{1}{2}\left\{\sigma_i^{%2b}\sigma_i^-,\rho\right\}\right)">

and the last term represent dephasing in <img src="https://render.githubusercontent.com/render/math?math=xy" style="vertical-align:middle"> plane:

<img src="https://render.githubusercontent.com/render/math?math=\mathcal{D}_2[\rho] = \sum_i g_{2,i} \left(\sigma_i^z \rho\sigma_i^z - \rho\right)">
