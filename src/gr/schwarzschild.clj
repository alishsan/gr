(ns gr.schwarzschild
  "Schwarzschild metric and geometry for General Relativity calculations."
  (:require [clojure.core.matrix :as m]))

(m/set-current-implementation :persistent-vector)

(defn schwarzschild-metric
  "Returns the Schwarzschild metric tensor g_μν at coordinates (t, r, θ, φ).
  
  Parameters:
  - r: radial coordinate
  - M: mass parameter (Schwarzschild radius = 2M)
  
  Returns a 4x4 matrix representing the metric tensor in (t, r, θ, φ) coordinates.
  The metric signature is (-, +, +, +)."
  [r M]
  (let [rs (* 2 M)  ; Schwarzschild radius
        f (- 1 (/ rs r))]
    (m/matrix [[(- f) 0 0 0]
               [0 (/ 1 f) 0 0]
               [0 0 (* r r) 0]
               [0 0 0 (* r r)]])))

(defn schwarzschild-metric-inverse
  "Returns the inverse Schwarzschild metric tensor g^μν."
  [r M]
  (let [rs (* 2 M)
        f (- 1 (/ rs r))]
    (m/matrix [[(/ -1 f) 0 0 0]
               [0 f 0 0]
               [0 0 (/ 1 (* r r)) 0]
               [0 0 0 (/ 1 (* r r))]])))

(defn vierbein
  "Computes the vierbein (tetrad) field e^a_μ for Schwarzschild spacetime.
  
  The vierbein relates the coordinate basis to an orthonormal basis.
  Returns a 4x4 matrix where rows are orthonormal basis indices (a) and
  columns are coordinate indices (μ).
  
  Assumes equatorial plane (θ = π/2) for simplicity."
  [r M]
  (let [rs (* 2 M)
        f (Math/sqrt (- 1 (/ rs r)))
        f-inv (Math/sqrt (/ 1 (- 1 (/ rs r))))]
    (m/matrix [[f 0 0 0]
               [0 f-inv 0 0]
               [0 0 r 0]
               [0 0 0 r]])))

(defn vierbein-inverse
  "Computes the inverse vierbein e^μ_a.
  
  Assumes equatorial plane (θ = π/2) for simplicity."
  [r M]
  (let [rs (* 2 M)
        f (Math/sqrt (- 1 (/ rs r)))
        f-inv (Math/sqrt (/ 1 (- 1 (/ rs r))))]
    (m/matrix [[(/ 1 f) 0 0 0]
               [0 f 0 0]
               [0 0 (/ 1 r) 0]
               [0 0 0 (/ 1 r)]])))

(defn christoffel-symbols
  "Computes Christoffel symbols Γ^λ_μν for Schwarzschild metric.
  
  Returns a 4x4x4 array where [λ μ ν] = Γ^λ_μν.
  Only non-zero components are computed."
  [r M]
  (let [rs (* 2 M)
        f (- 1 (/ rs r))
        df-dr (/ rs (* r r))
        r-sin-theta r]  ; Assuming θ = π/2 for equatorial plane
    (into-array
     (for [lambda (range 4)]
       (into-array
        (for [mu (range 4)]
          (into-array
           (for [nu (range 4)]
             (cond
               ;; Γ^t_tr = Γ^t_rt = (1/2) * (1/f) * df/dr
               (and (= lambda 0) (= mu 0) (= nu 1)) (* 0.5 (/ 1 f) df-dr)
               (and (= lambda 0) (= mu 1) (= nu 0)) (* 0.5 (/ 1 f) df-dr)
               ;; Γ^r_tt = (1/2) * f * df/dr
               (and (= lambda 1) (= mu 0) (= nu 0)) (* 0.5 f df-dr)
               ;; Γ^r_rr = -(1/2) * (1/f) * df/dr
               (and (= lambda 1) (= mu 1) (= nu 1)) (* -0.5 (/ 1 f) df-dr)
               ;; Γ^r_θθ = -r * f
               (and (= lambda 1) (= mu 2) (= nu 2)) (* -1 r f)
               ;; Γ^r_φφ = -r * f * sin²θ
               (and (= lambda 1) (= mu 3) (= nu 3)) (* -1 r f)
               ;; Γ^θ_rθ = Γ^θ_θr = 1/r
               (and (= lambda 2) (= mu 1) (= nu 2)) (/ 1 r)
               (and (= lambda 2) (= mu 2) (= nu 1)) (/ 1 r)
               ;; Γ^θ_φφ = -sinθ cosθ (assume θ = π/2, so = 0)
               (and (= lambda 2) (= mu 3) (= nu 3)) 0.0
               ;; Γ^φ_rφ = Γ^φ_φr = 1/r
               (and (= lambda 3) (= mu 1) (= nu 3)) (/ 1 r)
               (and (= lambda 3) (= mu 3) (= nu 1)) (/ 1 r)
               ;; Γ^φ_θφ = Γ^φ_φθ = cotθ (assume θ = π/2, so = 0)
               (and (= lambda 3) (= mu 2) (= nu 3)) 0.0
               (and (= lambda 3) (= mu 3) (= nu 2)) 0.0
               :else 0.0)))))))))

(defn spin-connection
  "Computes spin connection components ω^a_bμ for Schwarzschild spacetime.
  
  The spin connection is computed from the vierbein and Christoffel symbols.
  Returns a 4x4x4 array where [a b mu] = ω^a_bμ."
  [r M]
  (let [e (vierbein r M)
        e-inv (vierbein-inverse r M)
        gamma (christoffel-symbols r M)]
    (into-array
     (for [a (range 4)]
       (into-array
        (for [b (range 4)]
          (into-array
           (for [mu (range 4)]
             (reduce +
                     (for [nu (range 4)
                           lambda (range 4)]
                       (* (m/mget e a nu)
                          (m/mget e-inv lambda b)
                          (aget gamma lambda mu nu))))))))))))

