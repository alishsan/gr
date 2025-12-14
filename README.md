# gr - Dirac Equation in General Relativity

A Clojure library for solving the Dirac equation in curved spacetime, starting with the Schwarzschild metric.

## Overview

This project implements numerical solutions to the Dirac equation for spin-1/2 particles in General Relativity. The Dirac equation in curved spacetime takes the form:

```
(iγ^μ(∂_μ + Γ_μ) - m)ψ = 0
```

where:
- `γ^μ` are the curved-space Dirac matrices
- `Γ_μ` is the spin connection
- `ψ` is the Dirac spinor (4-component wave function)
- `m` is the particle mass

## Features

- **Schwarzschild Metric**: Complete implementation of the Schwarzschild metric and its geometric properties
- **Vierbein (Tetrad) Fields**: Orthonormal basis fields for converting between coordinate and local frames
- **Spin Connection**: Computation of spin connection components for the Dirac equation
- **Complex Number Support**: Full complex number arithmetic for proper treatment of the Dirac equation
- **Dirac Equation Solver**: Numerical integration methods (Euler and Runge-Kutta 4) for solving the Dirac equation
- **Dirac Current**: Computation of conserved currents with proper complex conjugation

## Installation

Add the following to your `project.clj` dependencies:

```clojure
[gr "0.1.0-SNAPSHOT"]
```

Or clone this repository and use it as a local dependency.

## Usage

### Basic Example: Schwarzschild Metric

```clojure
(require '[gr.schwarzschild :as sch])

;; Compute metric at radius r = 10M (where M is the Schwarzschild mass)
(def g (sch/schwarzschild-metric 10.0 1.0))
```

### Example: Solving the Dirac Equation

```clojure
(require '[gr.dirac :as dirac]
         '[gr.complex :as c :refer [complex-from-cartesian]]
         '[clojure.core.matrix :as m])

;; Initial conditions
(def r0 10.0)        ; Initial radius (must be > 2M)
(def r-final 20.0)   ; Final radius
(def dr 0.1)         ; Step size
(def M 1.0)          ; Schwarzschild mass parameter
(def m 0.1)          ; Particle mass

;; Initial spinor (can use real numbers or complex numbers)
(def psi0 [1.0 0.0 0.0 0.0])  ; Real spinor
;; Or with complex numbers:
;; (def psi0 [(complex-from-cartesian 1.0 0.0)
;;            (complex-from-cartesian 0.0 0.0)
;;            (complex-from-cartesian 0.0 0.0)
;;            (complex-from-cartesian 0.0 0.0)])

;; Solve using Runge-Kutta 4 method
(def solution (dirac/solve-dirac-schwarzschild r0 r-final dr psi0 M m :method :rk4))

;; Each element of solution is [r, psi] where psi is the spinor at radius r
```

### Example: Computing Dirac Current

```clojure
(require '[gr.dirac :as dirac]
         '[gr.complex :as c :refer [complex-from-cartesian re im]])

(def psi [1.0 0.0 0.0 0.0])  ; Spinor (real or complex)
(def r 10.0)                  ; Radial coordinate
(def M 1.0)                   ; Mass parameter

;; Compute conserved current j^μ
(def current (dirac/dirac-current psi r M))
;; Returns [j^t, j^r, j^θ, j^φ] (complex numbers)

;; Access real and imaginary parts
(println "j^t =" (re (nth current 0)) "+ i" (im (nth current 0)))
```

## Mathematical Background

### Schwarzschild Metric

The Schwarzschild metric in coordinates (t, r, θ, φ) is:

```
ds² = -(1 - 2M/r)dt² + (1 - 2M/r)⁻¹dr² + r²(dθ² + sin²θ dφ²)
```

where `M` is the mass parameter and `r_s = 2M` is the Schwarzschild radius.

### Vierbein Fields

The vierbein (tetrad) field `e^a_μ` relates the coordinate basis to an orthonormal basis:

```
g_μν = η_ab e^a_μ e^b_ν
```

where `η_ab` is the Minkowski metric.

### Spin Connection

The spin connection `ω^a_bμ` is computed from the vierbein and Christoffel symbols:

```
ω^a_bμ = e^a_ν (∂_μ e^ν_b + Γ^ν_μλ e^λ_b)
```

### Dirac Equation

The Dirac equation in curved spacetime requires:
1. Curved-space gamma matrices: `γ^μ = e^μ_a γ^a`
2. Spin connection term: `Γ_μ = (1/4) ω_abμ γ^a γ^b`
3. The full equation: `(iγ^μ(∂_μ + Γ_μ) - m)ψ = 0`

## Project Structure

```
gr/
├── src/
│   └── gr/
│       ├── core.clj          ; Main namespace with examples
│       ├── schwarzschild.clj ; Schwarzschild metric and geometry
│       └── dirac.clj         ; Dirac equation implementation
├── test/
│   └── gr/
│       └── core_test.clj     ; Unit tests
└── doc/
    └── intro.md              ; Introduction documentation
```

## Dependencies

- `org.clojure/clojure` - Clojure language
- `net.mikera/core.matrix` - Matrix operations
- `net.mikera/vectorz-clj` - Vector operations

## Complex Number Support

The library includes a full complex number implementation (`gr.complex`) with:
- Cartesian and polar coordinate representations
- Basic arithmetic operations (add, subtract, multiply, divide)
- Complex exponential and power functions
- Complex conjugation
- Complex integration functions

**Note on Matrix Operations**: The core matrix operations use real numbers for compatibility with `core.matrix`. Complex numbers are available for spinor components and can be used with the complex number library functions. For full complex matrix operations, users can work with real and imaginary parts separately or extend the matrix operations as needed.

## Limitations and Future Work

- Assumes equatorial plane (θ = π/2) for some calculations
- Numerical integration methods can be extended with adaptive step sizes
- Additional metrics (Kerr, Reissner-Nordström) can be added
- Visualization tools for solutions
- Enhanced matrix operations with complex numbers

## References

1. Dirac, P. A. M. (1928). "The quantum theory of the electron". Proceedings of the Royal Society of London.
2. Chandrasekhar, S. (1983). "The Mathematical Theory of Black Holes". Oxford University Press.
3. Parker, L., & Toms, D. J. (2009). "Quantum Field Theory in Curved Spacetime". Cambridge University Press.

## License

Copyright © 2025

This program and the accompanying materials are made available under the
terms of the Eclipse Public License 2.0 which is available at
http://www.eclipse.org/legal/epl-2.0.

This Source Code may also be made available under the following Secondary
Licenses when the conditions for such availability set forth in the Eclipse
Public License, v. 2.0 are satisfied: GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or (at your
option) any later version, with the GNU Classpath Exception which is available
at https://www.gnu.org/software/classpath/license.html.
