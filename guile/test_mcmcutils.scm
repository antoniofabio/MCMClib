(load "main.scm")

(use-modules (oop goops)
             (oop goops describe))

;;check gsl support functions
(gv2v (v2gv #(1 2 3)))
;#(1.0 2.0 3.0)
(gM2M (M2gM #2((1 2 3) (4 5 6))))
;#2((1.0 2.0 3.0) (4.0 5.0 6.0))

(define s (make-amh 'gauss-am
                    (new-gsl-rng (gsl-rng-default))
                    (make-guile-distrfun (lambda (x) 0.0))
                    (new-gsl-vector 1)
                    (M2gM #2((1.0)))
                    10))
;(do-ec (: i 100000) (update s)) ;;this should stop on values divergence

(define s (make-mh 'gauss-mrw
                    (new-gsl-rng (gsl-rng-default))
                    (make-guile-distrfun (lambda (x) 0.0))
                    (new-gsl-vector 1)
                    (M2gM #2((1.0)))))
(do-ec (: i 100000) (update s)) ;;this should not stop
(gv2v (slot-ref s 'x))
