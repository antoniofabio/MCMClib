(set! %load-path (cons "." %load-path))

(use-modules
 (oop goops)
 (oop goops describe)
 (ice-9 pretty-print)
 (srfi srfi-1)
 (srfi srfi-42)
 (swig gsl)
 (swig mcmclib)
 (swig mcmcutils))

(load "simul_aux.scm") ;auxiliary functions

(define dim 2) ;state space dimension
(define K 2) ;number of components of the fitted mixture

(define beta #(0.4 0.1 0.5))
(define mu (vector (make-filled-vector 0.0 dim)
                   (make-filled-vector 3.0 dim)
                   (make-filled-vector 8.0 dim)))
(define Sigma (vector-ec (: i 3) (diag 1.0 dim)))
(define pi-obj (make-guile-distrfun (make-mixnorm beta mu Sigma)))

(define x (new-gsl-vector dim))

(define (try-alpha alpha N)
  (let*
      ((trash (gsl-vector-set-all x 4.0))
       (s (make-amh 'raptor
                    (new-gsl-rng (gsl-rng-default)) ;random number generator
                    pi-obj ;target distribution object
                    x ;current chain value
                    5000 ;burnin
                    (diag 100.0 dim) ;starting proposal var-cov. matrix
                    (v2gv #(0.5 0.5)) ;beta-hat
                    (va2ca (vector (make-filled-vector -2.0 dim)
                                   (make-filled-vector 2.0 dim))) ;mu-hat
                    (ma2ca (vector-ec (: i K) (diag 1.0 dim))))) ;Sigma-hat
       (monitor-base '())
       (monitor-region '()))
    (mcmclib-raptor-set-alpha (slot-ref s 'c-ref) alpha)
    (do-ec (: i 5000) (update s)) ;burnin
    (set! monitor-base (new-mcmclib-monitor x))
    (set! monitor-region (make-monitor-region (vector->list mu) x))
    (do-ec (: i N)
           (begin
             (update s)
             (mcmclib-monitor-update monitor-base)
             (update monitor-region)))
    (mcmclib-monitor-fprintf-all monitor-base (stdout))
    (describe (slot-ref monitor-region 'trans-mat))))

(try-alpha 0.0 100000)
(try-alpha 0.5 100000)
(try-alpha 1.0 100000)
