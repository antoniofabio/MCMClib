(load "main.scm")

(define v (new-gsl-vector 1))
(gsl-vector-set v 0 3.0)
(test-distrfun (test-distrfun-alloc 5.0) v)
(test-distrfun (test-distrfun-alloc 2.0) v)

(define rng (new-gsl-rng (gsl-rng-default)))
(define lpdf-data (test-distrfun-alloc 1.0))
(define S (new-gsl-matrix 1 1))
(gsl-matrix-set-all S 0.1)
(gsl-vector-set v 0 0.5)
(define mh (mcmclib-gauss-mrw-alloc
            rng
            (test-distrfun-cb)
            lpdf-data
            v
            S))
(do-ec (: i 1e5)
       (mcmclib-mh-update mh))
(dv v)

(define mon (new-mcmclib-monitor v))
(define (update N)
  (do-ec (: i N)
         (begin
           (mcmclib-mh-update mh)
           (mcmclib-monitor-update mon))))

(update 1e5)
(mcmclib-monitor-fprintf-means mon (stdout))
(mcmclib-monitor-fprintf-all mon (stdout))

(define (f x)
  (if (< (gsl-vector-get x 0) 0)
      (log 0.0)
      0.0))

(define pf (guile-to-voidptr f))
(mcmclib-guile-lpdf pf v)

(define (dunif x)
  (let
      ((x0 (gsl-vector-get x 0)))
    (if
     (and
      (>= x0 0)
      (<= x0 1))
     0.0
     (log 0.0))))
(define mh (mcmclib-gauss-mrw-alloc
            rng
            (mcmclib-guile-lpdf-cb)
            (guile-to-voidptr dunif)
            v
            S))
(define mon (new-mcmclib-monitor v))
(update 10000)
(mcmclib-monitor-fprintf-all mon (stdout))