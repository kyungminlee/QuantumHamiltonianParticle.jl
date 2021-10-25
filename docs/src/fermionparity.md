# Notes on Fermion Parity

`QuantumHamiltonianParticle` uses Wigner-Jordan (WJ) transformatoin to represent fermion creation and annihilation operators
math```
c_{i}^{\dagger} = \sigma^{-} \prod_{j < i} \sigma^{z}_{i}
```
with the WJ string between site i and site 1.

The basis states are also defined using the σ⁻.
For example, the basis state |00110⟩ in a 5-site fermion system, with sites numbered 1 to 4 from the right, is defined as
math```
\vert 00110 \rangle \equiv \sigma_{2}^{-} \sigma_{3}^{-} \vert 00000 \rangle
```
Since the σ's on different sites commute with each other, the order of the σ⁻'s in the product does not matter.
This defines a natural convention for the states in terms of the *fermion* creation operators.
Since the WJ string runs toward site 1, applying creation operators to the vacuum state |00...0⟩ starting with larger site indices gives basis states consistent with the basis states
defined in terms of the σ⁻'s.
The basis state in the above examples can be written as
math```
\vert 00110 \rangle \equiv c_{2}^{\dagger} c_{3}^{\dagger} \vert 00000 \rangle.
```