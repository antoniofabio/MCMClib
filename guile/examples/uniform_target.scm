(use-modules
 (oop goops)
 (oop goops describe)
 (srfi srfi-1)
 (srfi srfi-42)
 (swig gsl)
 (swig mcmclib)
 (swig mcmcutils))

;;component-wise random walk on a (bounded) uniform distribution
(define (pi x)
  (let
      ((ans 0.0))
    (do-ec (: i (gsl-vector-size-get x))
           (if
            (or
             (< (gsl-vector-get x i) 0)
             (> (gsl-vector-get x i) 1))
            (set! ans (log 0.0))))
  ans))
(define rng (new-gsl-rng (gsl-rng-default)))
(define s (make-mh 'gauss-rw
                   rng
                   pi
                   (new-gsl-vector 5)
                   0.1))
(do-ec (: i 100000) (update s))

(define s (make-amh 'gauss-am
                   rng
                   pi
                   (new-gsl-vector 3)
                   (M2gM #(#(1 0 0) #(0 1 0) #(0 0 1)))
                   10))
(do-ec (: i 100000) (update s))
(gv2v (slot-ref s 'x))
