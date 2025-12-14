(ns gr.core-test
  (:require [clojure.test :refer :all]
            [gr.core :refer :all]
            [gr.schwarzschild :as sch]
            [gr.dirac :as dirac]
            [clojure.core.matrix :as m]))

(m/set-current-implementation :vectorz)

(deftest schwarzschild-metric-test
  (testing "Schwarzschild metric at large radius approaches Minkowski"
    (let [r 1000.0
          M 1.0
          g (sch/schwarzschild-metric r M)]
      ;; At large r, g_tt should approach -1, g_rr should approach 1
      (is (< (Math/abs (- (m/mget g 0 0) -1.0)) 0.01))
      (is (< (Math/abs (- (m/mget g 1 1) 1.0)) 0.01)))))

(deftest vierbein-test
  (testing "Vierbein and inverse vierbein relationship"
    (let [r 10.0
          M 1.0
          e (sch/vierbein r M)
          e-inv (sch/vierbein-inverse r M)
          product (m/mmul e e-inv)]
      ;; e * e-inv should be approximately identity
      (is (< (Math/abs (- (m/mget product 0 0) 1.0)) 0.001))
      (is (< (Math/abs (- (m/mget product 1 1) 1.0)) 0.001))
      (is (< (Math/abs (- (m/mget product 2 2) 1.0)) 0.001))
      (is (< (Math/abs (- (m/mget product 3 3) 1.0)) 0.001)))))

(deftest dirac-gamma-matrices-test
  (testing "Curved gamma matrices are computed correctly"
    (let [r 10.0
          M 1.0
          gamma-curved (dirac/curved-gamma-matrices r M)]
      ;; Should have 4 gamma matrices
      (is (= (count gamma-curved) 4))
      ;; Each should be a 4x4 matrix
      (is (= (m/shape (first gamma-curved)) [4 4])))))

(deftest dirac-solver-test
  (testing "Dirac equation solver produces results"
    (let [r0 10.0
          r-final 12.0
          dr 0.5
          M 1.0
          m 0.1
          psi0 [1.0 0.0 0.0 0.0]
          solution (dirac/solve-dirac-schwarzschild r0 r-final dr psi0 M m)]
      ;; Should produce multiple solution points
      (is (> (count solution) 1))
      ;; Each solution should be [r, psi] where psi is 4-component
      (let [[r psi] (first solution)
            psi-vec (vec (m/to-vector psi))]
        (is (number? r))
        (is (= (count psi-vec) 4))))))
