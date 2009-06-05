(load "main.scm")

(define v (new-gsl-vector 1))
(gsl-vector-set v 0 3.0)
(test-distrfun (test-distrfun-alloc 5.0) v)
(test-distrfun (test-distrfun-alloc 2.0) v)
