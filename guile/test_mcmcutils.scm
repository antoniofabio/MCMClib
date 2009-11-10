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

;;check gsl support functions
(gv2v (v2gv #(1 2 3)))
;#(1.0 2.0 3.0)
(gM2M (M2gM #(#(1 2 3) #(4 5 6))))
;#(#(1.0 2.0 3.0) #(4.0 5.0 6.0))

(define s (make-amh 'gauss-am
                    (new-gsl-rng (gsl-rng-default))
                    (make-guile-distrfun (lambda (x) 0.0))
                    (new-gsl-vector 1)
                    (M2gM #(#(1.0)))
                    10))
(do-ec (: i 100000) (update s)) ;;this should stop on values divergence

(define s (make-mh 'gauss-mrw
                    (new-gsl-rng (gsl-rng-default))
                    (make-guile-distrfun (lambda (x) 0.0))
                    (new-gsl-vector 1)
                    (M2gM #(#(1.0)))))
(do-ec (: i 100000) (update s)) ;;this should not stop
(gv2v (slot-ref s 'x))