;(chdir "/home/antonio/MCMClib/guile")
(set! %load-path (cons "." %load-path))

(use-modules
 (oop goops)
 (oop goops describe)
 (srfi srfi-1)
 (srfi srfi-42)
 (swig gsl)
 (swig mcmclib)
 (swig mcmcutils))

(define s (make-amh 'gauss-am
                    (new-gsl-rng (gsl-rng-default))
                    (make-guile-distrfun (lambda (x) 0.0))
                    (new-gsl-vector 1)
                    (let ((S (new-gsl-matrix 1 1))) (gsl-matrix-set S 0 0 1.0) S)
                    10))
(do-ec (: i 100000) (update s)) ;;this should stop on values divergence

(define s (make-mh 'gauss-mrw
                    (new-gsl-rng (gsl-rng-default))
                    (make-guile-distrfun (lambda (x) 0.0))
                    (new-gsl-vector 1)
                    (let ((S (new-gsl-matrix 1 1))) (gsl-matrix-set S 0 0 1.0) S)))
(do-ec (: i 100000) (update s)) ;;this should not stop
(gv2v (slot-ref s 'x))
