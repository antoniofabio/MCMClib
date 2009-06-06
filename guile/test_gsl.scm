(set! %load-path (cons "." %load-path))

(use-modules (srfi srfi-1)
	     (srfi srfi-42)
             (swig gsl))

(define v (new-gsl-vector 3))
(do-ec (: i 3)
       (gsl-vector-set v i (* i i)))
(gsl-vector-fprintf (stdout) v "%f")

(define fa (new-vectorArray 3))
(define (va2ca va)
  (let*
      ((n (vector-length va))
       (ca (new-vectorArray n)))
    (do-ec (: i n)
           (vectorArray-setitem ca i (vector-ref va i)))
    ca))
(va2ca
 (vector
  (new-gsl-vector 2)
  (new-gsl-vector 3)))

(define m (new-gsl-matrix 3 3))
(do-ec (: i 3) (: j 3)
       (gsl-matrix-set m i j (* i j)))
(gsl-matrix-fprintf (stdout) m "%f")

(define fa (new-matrixArray 3))
(define (ma2ca ma)
  (let*
      ((n (vector-length ma))
       (ca (new-matrixArray n)))
    (do-ec (: i n)
           (matrixArray-setitem ca i (vector-ref ma i)))
    ca))
(ma2ca
 (vector
  (new-gsl-matrix 2 2)
  (new-gsl-matrix 3 3)))

(define rng (new-gsl-rng (gsl-rng-default)))
(define lr )
(define (mean l)
  (/
   (fold (lambda (x y) (+ x y)) 0 l)
   (length l)))
(mean (list-ec (: i 100) (gsl-rng-uniform rng)))
