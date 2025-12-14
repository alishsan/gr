# Introduction to gr

## What is gr?

`gr` is a Clojure library for solving the Dirac equation in curved spacetime, with initial support for the Schwarzschild metric. The Dirac equation describes the quantum mechanical behavior of spin-1/2 particles (like electrons) in the presence of gravitational fields.

## Why Solve the Dirac Equation in GR?

The Dirac equation in curved spacetime is fundamental to understanding:

- **Quantum field theory in curved spacetime**: How quantum fields behave near black holes
- **Hawking radiation**: The theoretical prediction that black holes emit radiation
- **Particle physics in strong gravitational fields**: Behavior of fermions near compact objects
- **Unruh effect**: The relationship between acceleration and temperature in quantum field theory

## The Schwarzschild Metric

The Schwarzschild metric describes the spacetime geometry around a non-rotating, uncharged black hole. It's the simplest exact solution to Einstein's field equations and serves as an excellent starting point for understanding curved spacetime effects.

Key features:
- **Event horizon**: At `r = 2M` (Schwarzschild radius)
- **Singularity**: At `r = 0`
- **Asymptotically flat**: Approaches Minkowski spacetime at large distances

## How It Works

### 1. Metric and Geometry

The library computes the Schwarzschild metric tensor `g_μν` and related geometric quantities:

```clojure
(require '[gr.schwarzschild :as sch])

;; Metric at radius r = 10M
(def g (sch/schwarzschild-metric 10.0 1.0))
```

### 2. Vierbein Fields

To work with spinors in curved spacetime, we need vierbein (tetrad) fields that provide an orthonormal basis:

```clojure
(def e (sch/vierbein 10.0 1.0))
```

### 3. Spin Connection

The spin connection encodes how spinors are parallel transported in curved spacetime:

```clojure
(def omega (sch/spin-connection 10.0 1.0))
```

### 4. Dirac Equation

The full Dirac equation combines:
- Curved-space gamma matrices
- Spin connection terms
- The Dirac spinor field

```clojure
(require '[gr.dirac :as dirac])

;; Solve the equation numerically
(def solution (dirac/solve-dirac-schwarzschild 
                10.0 20.0 0.1  ; r0, r-final, dr
                [1.0 0.0 0.0 0.0]  ; initial spinor
                1.0 0.1  ; M, m
                :method :rk4))
```

## Physical Interpretation

The Dirac spinor `ψ` has 4 components representing different spin and energy states. The conserved current `j^μ = ψ̄ γ^μ ψ` represents the probability current density, which must be conserved (∇_μ j^μ = 0) in the absence of sources.

## Getting Started

1. **Install dependencies**: `lein deps`
2. **Start a REPL**: `lein repl`
3. **Try the examples**: See the README for usage examples
4. **Run tests**: `lein test`

## Next Steps

- Explore different initial conditions
- Study behavior near the event horizon
- Compute scattering cross-sections
- Extend to other metrics (Kerr, Reissner-Nordström)

## Further Reading

- "Quantum Field Theory in Curved Spacetime" by Parker & Toms
- "The Mathematical Theory of Black Holes" by Chandrasekhar
- "Quantum Fields in Curved Space" by Birrell & Davies
