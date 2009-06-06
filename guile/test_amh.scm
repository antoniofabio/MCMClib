(load "main.scm")

(define rng (new-gsl-rng (gsl-rng-default)))
(define v (new-gsl-vector 1))
(gsl-vector-set v 0 0.1)
(define rng (new-gsl-rng (gsl-rng-default)))
(define S (new-gsl-matrix 1 1))
(gsl-matrix-set-all S 0.1)
(define (dunif x)
  (let
      ((x0 (gsl-vector-get x 0)))
    (if
     (and
      (>= x0 0)
      (<= x0 1))
     0.0
     (log 0.0))))
(define amh (mcmclib-gauss-am-alloc
             rng
             (mcmclib-guile-lpdf-cb)
             (guile-to-voidptr dunif)
             v
             S
             100))

(define mon (new-mcmclib-monitor v))
(define (update N)
  (do-ec (: i N)
         (begin
           (mcmclib-amh-update amh)
           (mcmclib-monitor-update mon))))
(update 10000)
(mcmclib-monitor-fprintf-all mon (stdout))
