(set! %load-path (cons "." %load-path))

(use-modules (srfi srfi-1)
	     (srfi srfi-42)
             (swig gsl))

(define v (new-gsl-vector 3))
(do-ec (: i 3)
       (gsl-vector-set v i (* i i)))
(gsl-vector-fprintf (stdout) v "%f")

(define m (new-gsl-matrix 3 3))
(do-ec (: i 3) (: j 3)
       (gsl-matrix-set m i j (* i j)))
(gsl-matrix-fprintf (stdout) m "%f")

(define rng (new-gsl-rng (gsl-rng-default)))
(define lr )
(define (mean l)
  (/
   (fold (lambda (x y) (+ x y)) 0 l)
   (length l)))
(mean (list-ec (: i 100) (gsl-rng-uniform rng)))
