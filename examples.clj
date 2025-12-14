(ns examples
  "Example usage of the gr library for solving Dirac equation in Schwarzschild spacetime."
  (:require [gr.schwarzschild :as sch]
            [gr.dirac :as dirac]
            [gr.complex :as c :refer [complex-from-cartesian re im]]
            [clojure.core.matrix :as m]))

(m/set-current-implementation :vectorz)

;; Example 1: Compute Schwarzschild metric
(println "Example 1: Schwarzschild Metric")
(println "================================")
(let [r 10.0  ; 10 times the Schwarzschild radius
      M 1.0   ; Mass parameter
      g (sch/schwarzschild-metric r M)]
  (println "Metric at r = 10M:")
  (println (m/pm g))
  (println))

;; Example 2: Vierbein fields
(println "Example 2: Vierbein Fields")
(println "===========================")
(let [r 10.0
      M 1.0
      e (sch/vierbein r M)
      e-inv (sch/vierbein-inverse r M)]
  (println "Vierbein e^a_μ:")
  (println (m/pm e))
  (println "Inverse vierbein e^μ_a:")
  (println (m/pm e-inv))
  (println "Check: e * e-inv should be identity")
  (println (m/pm (m/mmul e e-inv)))
  (println))

;; Example 3: Curved-space gamma matrices
(println "Example 3: Curved-Space Gamma Matrices")
(println "=======================================")
(let [r 10.0
      M 1.0
      gamma-curved (dirac/curved-gamma-matrices r M)]
  (println "γ^t (time component):")
  (println (m/pm (nth gamma-curved 0)))
  (println "γ^r (radial component):")
  (println (m/pm (nth gamma-curved 1)))
  (println))

;; Example 4: Solve Dirac equation with complex numbers
(println "Example 4: Solving Dirac Equation")
(println "===================================")
(let [r0 10.0        ; Start at 10M
      r-final 15.0   ; End at 15M
      dr 0.5         ; Step size
      M 1.0          ; Schwarzschild mass
      m 0.1          ; Particle mass
      ;; Initial spinor with complex numbers
      psi0 [(complex-from-cartesian 1.0 0.0)
            (complex-from-cartesian 0.0 0.0)
            (complex-from-cartesian 0.0 0.0)
            (complex-from-cartesian 0.0 0.0)]
      solution (dirac/solve-dirac-schwarzschild r0 r-final dr psi0 M m :method :rk4)]
  (println "Solution points:" (count solution))
  (println "First few points:")
  (doseq [[r psi] (take 3 solution)]
    (let [psi-vec (m/to-vector psi)]
      (println (format "r = %.2f" r))
      (doseq [[i comp] (map-indexed vector psi-vec)]
        (println (format "  psi[%d] = %.4f + %.4fi" i (re comp) (im comp)))))
    (println))
  (println))

;; Example 5: Dirac current with complex numbers
(println "Example 5: Dirac Current")
(println "=========================")
(let [psi [(complex-from-cartesian 1.0 0.0)
           (complex-from-cartesian 0.0 0.0)
           (complex-from-cartesian 0.0 0.0)
           (complex-from-cartesian 0.0 0.0)]
      r 10.0
      M 1.0
      current (dirac/dirac-current psi r M)]
  (println "Current j^μ (complex numbers):")
  (doseq [[mu j] (map-indexed vector current)]
    (println (format "j^%d = %.4f + %.4fi" mu (re j) (im j))))
  (println))

