(ns gr.core
  "Main namespace for General Relativity calculations.
  
  This library provides tools for solving the Dirac equation in curved spacetime,
  starting with the Schwarzschild metric."
  (:require [gr.schwarzschild :as sch]
            [gr.dirac :as dirac]
            [gr.complex :as c :refer [complex-from-cartesian]]))

(defn example-schwarzschild-metric
  "Example: Compute Schwarzschild metric at a given radius."
  [r M]
  (sch/schwarzschild-metric r M))

(defn example-dirac-solution
  "Example: Solve Dirac equation in Schwarzschild spacetime.
  
  Parameters:
  - r0: initial radius (must be > 2M to avoid singularity)
  - r-final: final radius
  - dr: step size
  - M: Schwarzschild mass
  - m: particle mass
  - psi0: initial spinor [a b c d] (can be real numbers or complex numbers)"
  [r0 r-final dr M m psi0]
  (dirac/solve-dirac-schwarzschild r0 r-final dr psi0 M m))

(defn example-complex-spinor
  "Example: Create a complex spinor for the Dirac equation.
  
  Returns a 4-component spinor with complex numbers."
  []
  [(complex-from-cartesian 1.0 0.0)  ; Real component
   (complex-from-cartesian 0.0 0.5)  ; Imaginary component
   (complex-from-cartesian 0.0 0.0)
   (complex-from-cartesian 0.0 0.0)])

(defn example-current
  "Example: Compute Dirac current for a given spinor."
  [psi r M]
  (dirac/dirac-current psi r M))
