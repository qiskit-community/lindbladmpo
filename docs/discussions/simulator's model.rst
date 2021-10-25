#################
Lindblad dynamics
#################

The state of an open quantum system is defined by a density matrix :math:`\rho`.
In the case where the system is made of of qubits (two level systems), the Hilbert space has dimension 2^N (therefore the density matrix is 2^N by 2^N in size).

When the system is coupled to the environment (memory-less bath), it's state evolves over time.
The notation of such evolution is :math:`\mathcal{L}[\rho]`.
This Liouvillian is a (linear) superoperator that maps the operator :math:`\rho` to another operator - the time derivative of :math:`\rho`, and is often described using the master equation:

.. math::

            \frac{\partial}{\partial t}\rho = \mathcal{L}[\rho]\equiv -\frac{i}{\hbar}[H,\rho]+\mathcal{D}[\rho]

where :math:`[\cdot,\cdot]` describes the commutator of two operators.

The first part of the equation, :math:`\frac{i}{\hbar}[H,\rho]`, expresses the unitary evolution due to interactions and coherent driving terms,
while :math:`\mathcal{D}[\rho]`  (sometimes known as the dissipator, or Lindbladian) is a superoperator that accounts for incoherent dephasing and relaxation processes due to the environment.

#############
The simulator
#############

The simulator simulates numerically Lindbladian dynamics and offers approximate solution to the master equation.
The simulator's solver is based on matrix product states (MPS) and matrix product operators (MPO).
The matrix product state is used to represent the state of the system, while the matrix product operator is used to represent the time-evolution operator.
The solver integrates the dynamics in fixed time steps :math:`\tau` defined by the user.

The Model
_________

The simulator can model interacting qubits in a planar lattice with an arbitrary connectivity defined by the user.

The simulator uses a fixed Hamiltonian with coefficients which can be defined by the user. The Hamiltonian is defined in the rotating frame with
respect to the frequency of an identical drive applied to all qubits (which are not necessarily identical themselves, however, and may be driven at different amplitudes).
The Hamiltonian is represented as the sum of on-site terms and an interaction part T:

.. math::

            \frac{H}{\hbar} = \sum_{i}\frac{1}{2}\left[h_{z,i}\sigma_i^z  + h_{x,i}\sigma_i^x + h_{y,i}\sigma_i^y\right] + T

with the interaction being:

.. math::
            T = \sum_{ i}^N\sum_{ j\neq i}^N \left(J_{ij}\sigma^+_i \sigma^-_{j} +{\rm h.c.} +\frac{1}{2} J_{ij}^z \sigma^z_i \sigma^z_{j}\right)=\\ \frac{1}{2}\sum_{ i}^N\sum_{ j\neq i}^N \left(J_{ij}\sigma^x_i \sigma^x_{j} + J_{ij}\sigma^y_i \sigma^y_{j} + J_{ij}^z \sigma^z_i \sigma^z_{j}\right)

while :math:`\sigma_i^a, a=\{x,y,z\}` are the Pauli matrices at each site and :math:`\sigma^{\pm}_i = {\sigma^{x}_i\pm i\sigma^y_i}/{2}` are the ladder operators.

The Lindbladian used in the simulator is also fixed with coefficients which can be defined by the user. The Lindbladian is a sum of three terms -

.. math::

            \mathcal{D} = \sum_j \mathcal{D}_j

The first two terms represent exchange with a thermal bath:

.. math::

            \mathcal{D}_0[\rho] = \sum_i g_{0,i}\left(\sigma_i^+ \rho\sigma_i^- - \frac{1}{2} \{\sigma_i^- \sigma_i^+,\rho\}\right)

            \mathcal{D}_1[\rho]=\sum_i g_{1,i}\left( \sigma_i^-\rho \sigma_i^{+}-\frac{1}{2}\left\{\sigma_i^{+}\sigma_i^-,\rho\right\}\right)

and the last term represent dephasing in :math:`xy` plane:

.. math::

            \mathcal{D}_2[\rho] = \sum_i g_{2,i} \left(\sigma_i^z \rho\sigma_i^z - \rho\right)