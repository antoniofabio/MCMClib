(use-modules (srfi srfi-42)
	     (swig gsl-utils))
(use-modules ((swig gsl) :renamer gsl-renamer))

(define v (sv->gv #(2 1 -2)))
(call-with-output-file "/dev/stdout"
  (lambda (port)
    (gv-fprintf port v "%.3f")))

(define m (make-gm-diag 3 2.5))
(gm->sm m)

(define rng (new-gsl-rng (gsl-rng-default)))
(gsl-rng-uniform rng)

(vector-ec (:gv xi v) xi)
(define v (sv->gv #(11 10 12)))
(do-ec (:gv xi (index i) v) (format #t "v[~a] = ~a\n" i xi))
(do-ec (:gv-along i v) (format #t "v[~a] = ~a\n" i (gv-get v i)))
