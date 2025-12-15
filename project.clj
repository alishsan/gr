(defproject gr "0.1.0-SNAPSHOT"
  :description "Dirac equation solver in General Relativity, starting with Schwarzschild spacetime"
  :url "https://github.com/alishsan/gr"
  :license {:name "EPL-2.0 OR GPL-2.0-or-later WITH Classpath-exception-2.0"
            :url "https://www.eclipse.org/legal/epl-2.0/"}
  :dependencies [[org.clojure/clojure "1.10.3"]
                 [net.mikera/core.matrix "0.63.0"]
                 [generateme/fastmath "3.0.0-alpha4-SNAPSHOT" :exclusions [com.github.haifengl/smile-mkl]]]
  :repl-options {:init-ns gr.core})
