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

;;Target distribution
(define beta #(0.4 0.1 0.5)) ;weights
(define mu (vector (make-filled-vector 0.0 dim)
                   (make-filled-vector 3.0 dim)
                   (make-filled-vector 8.0 dim))) ;means
(define Sigma (vector-ec (: i 3) (diag 1.0 dim))) ;variances
(define pi-obj (make-guile-distrfun (make-mixnorm beta mu Sigma))) ;distribution object

(define x (new-gsl-vector dim)) ;current chain value

(define (try-alpha alpha N)
  "run RAPTOR a simulation with alpha='alpha', 'N' iterations after burnin.
   At the end, prints some diagnostics on screen"
  (let*
      ((trash (gsl-vector-set-all x 4.0))
       (s (make-amh 'raptor
                    (new-gsl-rng (gsl-rng-default)) ;random number generator
                    pi-obj ;target distribution object
                    x ;current chain value
                    5000 ;burnin
                    (diag 100.0 dim) ;starting proposal var-cov. matrix
                    (v2gv #(0.5 0.5)) ;beta-hat
                    (va2ca (vector (make-filled-vector -10.0 dim)
                                   (make-filled-vector 10.0 dim))) ;mu-hat
                    (ma2ca (vector-ec (: i K) (diag 1.0 dim))))) ;Sigma-hat
       (monitor-base '())
       (monitor-region '()))
    (mcmclib-raptor-set-alpha (slot-ref s 'c-ref) alpha)
    ;(mcmclib-raptor-set-alpha-fun-identity (slot-ref s 'c-ref))
    (do-ec (: i 5000) (update s)) ;burnin
    (set! monitor-base (new-mcmclib-monitor x))
    (set! monitor-region (make-monitor-region (vector->list mu) x))
    (do-ec (: i N)
           (begin
             (update s)
             (mcmclib-monitor-update monitor-base)
             (update monitor-region)))
    (mcmclib-monitor-fprintf-all monitor-base (stdout))
    (describe (slot-ref monitor-region 'trans-mat))
    (pretty-print (gv2v (mcmclib-raptor-gamma-beta-hat-get
                         (mcmclib-raptor-gamma-get (slot-ref s 'c-ref)))))
    (do-ec (: k K) (pretty-print (gv2v (vectorArray-getitem (mcmclib-raptor-gamma-mu-hat-get
                                                             (mcmclib-raptor-gamma-get (slot-ref s 'c-ref))) k))))
    (do-ec (: k K) (pretty-print (gM2M (matrixArray-getitem (mcmclib-raptor-gamma-Sigma-hat-get
                                                             (mcmclib-raptor-gamma-get (slot-ref s 'c-ref))) k))))))

(try-alpha 0.0 10000)
(try-alpha 0.5 10000)
(try-alpha 1.0 10000)
