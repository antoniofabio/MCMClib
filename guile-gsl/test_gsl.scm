(use-modules (srfi srfi-42)
	     (swig gsl-utils))
(use-modules ((swig gsl) :renamer gsl-renamer))

(define v (new-gsl-vector 3))
(gsl-vector-set v 0 2.0)
(gsl-vector-get v 0)
(gsl-vector-set v 1 1.0)
(gsl-vector-set v 2 -2.0)
(call-with-output-file "/dev/stdout"
  (lambda (port)
    (gsl-vector-fprintf port v "%.3f")))

(define m (new-gsl-matrix 3 3))
(gsl-matrix-set m 0 0 2.0)
(gsl-matrix-get m 0 0)

(define rng (new-gsl-rng (gsl-rng-default)))
(gsl-rng-uniform rng)

(vector-ec (:gv xi v) xi)
(define v (sv->gv #(11 10 12)))
(do-ec (:gv xi (index i) v) (format #t "v[~a] = ~a\n" i xi))
(do-ec (:gv-along i v) (format #t "v[~a] = ~a\n" i (gsl-vector-get v i)))
