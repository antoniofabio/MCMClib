(load "main.scm")

(define rng (new-gsl-rng (gsl-rng-default)))
(define v (new-gsl-vector 1))
(gsl-vector-set v 0 0.1)
(define S (diag 1.0 1))
(define Sk (vector-ec (: k 2) S))
(define (reg x)
  (if (< (gsl-vector-get x 0) 0) 0 1))
(define (dunif x)
  (let
      ((x0 (gsl-vector-get x 0)))
    (if
     (and
      (>= x0 0)
      (<= x0 1))
     0.0
     (log 0.0))))
(define s (make-rapt
           rng
           dunif
           v
           100 ;t0
           S ;sigma-whole
           2 ;K
           Sk;sigma-local
           (make-guile-regionfun reg)))

(do-ec (: i 10000) (update s))
