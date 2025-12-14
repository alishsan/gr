# Change Log
All notable changes to this project will be documented in this file. This change log follows the conventions of [keepachangelog.com](http://keepachangelog.com/).

## [Unreleased]

## [0.1.0] - 2025-01-XX
### Added
- Initial implementation of Dirac equation solver in General Relativity
- Schwarzschild metric computation (`gr.schwarzschild`)
  - Metric tensor computation
  - Inverse metric tensor
  - Vierbein (tetrad) fields and their inverses
  - Christoffel symbols
  - Spin connection components
- Dirac equation implementation (`gr.dirac`)
  - Flat-space gamma matrices (Dirac representation)
  - Curved-space gamma matrices for Schwarzschild spacetime
  - Spin connection terms for Dirac equation
  - Numerical solver with Euler and Runge-Kutta 4 methods
  - Dirac current computation
- Comprehensive documentation
  - README with usage examples and mathematical background
  - Introduction guide in `doc/intro.md`
  - Example usage file
- Unit tests for core functionality

### Dependencies
- `org.clojure/clojure` "1.10.3"
- `net.mikera/core.matrix` "0.63.0"
- `net.mikera/vectorz-clj` "0.48.0"
