(ns gr.dirac
  "Dirac equation in curved spacetime.
  
  Implements the Dirac equation for spin-1/2 particles in General Relativity.
  The equation takes the form:
  
  (iγ^μ(∂_μ + Γ_μ) - m)ψ = 0
  
  where γ^μ are the curved-space Dirac matrices, Γ_μ is the spin connection,
  and ψ is the Dirac spinor."
  (:require [clojure.core.matrix :as m]
            [gr.schwarzschild :as sch]
            [gr.complex :as c :refer [complex-from-cartesian re im add mul subt]]))

(m/set-current-implementation :vectorz)

;; Flat-space Dirac matrices (gamma matrices) in Dirac representation
(def gamma-0-flat
  "γ^0 in flat spacetime (Dirac representation)"
  (m/matrix [[1 0 0 0]
             [0 1 0 0]
             [0 0 -1 0]
             [0 0 0 -1]]))

(def gamma-1-flat
  "γ^1 in flat spacetime"
  (m/matrix [[0 0 0 1]
             [0 0 1 0]
             [0 -1 0 0]
             [-1 0 0 0]]))

(def gamma-2-flat
  "γ^2 in flat spacetime.
  Note: In the full Dirac equation, this matrix is multiplied by i (the imaginary unit).
  For core.matrix compatibility, we store the real part only.
  γ^2 = i * [[0 0 0 -1]
             [0 0 1 0]
             [0 1 0 0]
             [-1 0 0 0]]"
  (m/matrix [[0 0 0 -1]
             [0 0 1 0]
             [0 1 0 0]
             [-1 0 0 0]]))

(def gamma-3-flat
  "γ^3 in flat spacetime"
  (m/matrix [[0 0 1 0]
             [0 0 0 -1]
             [-1 0 0 0]
             [0 1 0 0]]))

(def flat-gamma-matrices [gamma-0-flat gamma-1-flat gamma-2-flat gamma-3-flat])

;; Note: Complex number support is available via gr.complex for spinor components.
;; Matrix operations use real numbers for compatibility with core.matrix.
;; Users can work with complex spinors by handling real and imaginary parts separately.

(defn curved-gamma-matrices
  "Computes curved-space Dirac matrices γ^μ = e^μ_a γ^a.
  
  Parameters:
  - r: radial coordinate
  - M: mass parameter
  
  Returns a vector of 4 matrices [γ^t, γ^r, γ^θ, γ^φ]
  Note: γ^2 contains an implicit factor of i that must be handled in operations."
  [r M]
  (let [e-inv (sch/vierbein-inverse r M)]
    (mapv (fn [mu]
            (reduce m/add
                    (map (fn [a]
                           (let [coeff (m/mget e-inv mu a)
                                 gamma-a (nth flat-gamma-matrices a)]
                             (m/mul coeff gamma-a)))
                         (range 4))))
          (range 4))))

(defn spin-connection-term
  "Computes the spin connection term Γ_μ for the Dirac equation.
  
  Γ_μ = (1/4) * ω_abμ * γ^a * γ^b
  
  Returns a 4x4 matrix representing Γ_μ."
  [r M mu]
  (let [omega (sch/spin-connection r M)
        gamma-flat flat-gamma-matrices]
    (reduce m/add
            (for [a (range 4)
                  b (range 4)]
              (let [omega-ab-mu (aget omega a b mu)
                    gamma-a (nth gamma-flat a)
                    gamma-b (nth gamma-flat b)]
                (m/mul (* 0.25 omega-ab-mu)
                       (m/mmul gamma-a gamma-b)))))))

(defn dirac-equation-rhs
  "Computes the right-hand side of the Dirac equation for numerical integration.
  
  The Dirac equation: (iγ^μ(∂_μ + Γ_μ) - m)ψ = 0
  Rearranged for time evolution: ∂_t ψ = -i(γ^0)^(-1) * [γ^j(∂_j + Γ_j) - im]ψ
  
  Parameters:
  - psi: Dirac spinor (4-component complex vector)
  - r: radial coordinate
  - M: mass parameter
  - m: particle mass
  - dr: radial step size (for computing ∂_r numerically)
  
  Returns dψ/dt as a 4-component vector."
  [psi r M m dr]
  (let [gamma-curved (curved-gamma-matrices r M)
        gamma-0 (nth gamma-curved 0)
        gamma-r (nth gamma-curved 1)
        gamma-theta (nth gamma-curved 2)
        gamma-phi (nth gamma-curved 3)
        gamma-0-inv (m/inverse gamma-0)
        psi-vec (m/array psi)
        
        ;; Spin connection terms
        gamma-conn-0 (spin-connection-term r M 0)
        gamma-conn-r (spin-connection-term r M 1)
        gamma-conn-theta (spin-connection-term r M 2)
        gamma-conn-phi (spin-connection-term r M 3)
        
        ;; For simplicity, assume spherical symmetry and time-independence
        ;; ∂_r ψ approximated as zero for now (can be improved with finite differences)
        dpsi-dr (m/zero-vector 4)
        
        ;; Compute spatial derivative terms
        spatial-term (m/add
                      (m/mmul gamma-r (m/add dpsi-dr (m/mmul gamma-conn-r psi-vec)))
                      (m/mmul gamma-theta (m/mmul gamma-conn-theta psi-vec))
                      (m/mmul gamma-phi (m/mmul gamma-conn-phi psi-vec)))
        
        ;; Mass term: -im * psi (where i is the imaginary unit)
        ;; Note: Full complex number support would multiply by complex number (0, -m)
        ;; For now, we work with real matrices for core.matrix compatibility
        mass-term (m/mul m psi-vec)
        
        ;; Total spatial + mass contribution
        rhs-spatial (m/sub spatial-term mass-term)]
    
    ;; Multiply by -i * (γ^0)^(-1)
    ;; Note: The -i factor is handled implicitly in the simplified real version
    ;; Full implementation would multiply by complex number (0, -1)
    (m/mmul gamma-0-inv (m/mul -1 rhs-spatial))))

(defn solve-dirac-schwarzschild
  "Solves the Dirac equation in Schwarzschild spacetime using numerical integration.
  
  Parameters:
  - r0: initial radial coordinate
  - r-final: final radial coordinate
  - dr: radial step size
  - psi0: initial Dirac spinor (4-component vector)
  - M: Schwarzschild mass parameter
  - m: particle mass
  - method: integration method (:euler or :rk4, default :rk4)
  
  Returns a sequence of [r, psi] pairs."
  [r0 r-final dr psi0 M m & {:keys [method] :or {method :rk4}}]
  (let [rk4-step (fn [r psi dr]
                   (let [k1 (dirac-equation-rhs psi r M m dr)
                         k2 (dirac-equation-rhs (m/add psi (m/mul (* 0.5 dr) k1)) 
                                                (+ r (* 0.5 dr)) M m dr)
                         k3 (dirac-equation-rhs (m/add psi (m/mul (* 0.5 dr) k2))
                                                (+ r (* 0.5 dr)) M m dr)
                         k4 (dirac-equation-rhs (m/add psi (m/mul dr k3))
                                                (+ r dr) M m dr)]
                     (m/add psi (m/mul (/ dr 6.0) 
                                       (m/add k1 (m/mul 2 k2) (m/mul 2 k3) k4)))))
        
        euler-step (fn [r psi dr]
                     (m/add psi (m/mul dr (dirac-equation-rhs psi r M m dr))))
        
        step-fn (case method
                  :rk4 rk4-step
                  :euler euler-step
                  rk4-step)]
    
    (loop [r r0
           psi (m/array psi0)
           result []]
      (if (>= r r-final)
        (reverse (conj result [r psi]))
        (let [psi-new (step-fn r psi dr)]
          (recur (+ r dr) psi-new (conj result [r psi])))))))

(defn dirac-current
  "Computes the conserved Dirac current j^μ = ψ̄ γ^μ ψ.
  
  Parameters:
  - psi: Dirac spinor (vector, can contain real numbers or complex numbers)
  - r: radial coordinate
  - M: mass parameter
  
  Returns a 4-vector [j^t, j^r, j^θ, j^φ].
  
  ψ̄ = ψ† γ^0 where † denotes Hermitian conjugate (complex conjugate transpose).
  
  Note: For complex spinors, use gr.complex/complex-conjugate for proper conjugation.
  This version works with real spinors for core.matrix compatibility."
  [psi r M]
  (let [psi-vec (m/array psi)
        ;; For real spinors, transpose is sufficient
        ;; For complex spinors, would need complex conjugate transpose
        psi-transpose (m/transpose psi-vec)
        gamma-0 (nth flat-gamma-matrices 0)
        psi-bar (m/mmul psi-transpose gamma-0)
        gamma-curved (curved-gamma-matrices r M)]
    (mapv (fn [mu]
            (m/mmul psi-bar (m/mmul (nth gamma-curved mu) psi-vec)))
          (range 4))))

