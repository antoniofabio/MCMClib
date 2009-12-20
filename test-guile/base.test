(load "main.scm")

(define v (make-filled-vector 3.0 1))

(define rng (new-gsl-rng (gsl-rng-default)))
(define S (new-gsl-matrix 1 1))
(gsl-matrix-set-all S 0.1)
(gsl-vector-set v 0 0.5)
(define mh (mcmclib-gauss-mrw-alloc
            rng
            (lambda (x) (if (< (abs (gsl-vector-get x 0)) 1.0) 0.0 (log 0.0)))
            v
            S))
(do-ec (: i 1e5)
       (mcmclib-mh-update mh))
(gv2v v)

(define mon (new-mcmclib-monitor v))
(define (update N)
  (do-ec (: i N)
         (begin
           (mcmclib-mh-update mh)
           (mcmclib-monitor-update mon))))

(update 1e5)
(mcmclib-monitor-fprintf-means mon (stdout))
(mcmclib-monitor-fprintf-all mon (stdout))

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
            dunif
            v
            S))
(define mon (new-mcmclib-monitor v))
(update 10000)
(mcmclib-monitor-fprintf-all mon (stdout))
