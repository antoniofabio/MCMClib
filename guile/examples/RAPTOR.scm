(set! %load-path (cons "." %load-path))

(use-modules (srfi srfi-1)
	     (srfi srfi-42)
             (swig gsl)
             (swig mcmclib))

(define *n* 5) ;;target distribution dimension
(define *K* 2) ;;number of mixture components
(define *T0* 100) ;;burn-in length
(define *N* 1000) ;;number of iterations after burn-in
(define *rng* (new-gsl-rng (gsl-rng-default))) ;;random number generator

(define *beta-hat* (new-gsl-vector *K*)) ;;starting mixture weigths estimates
(gsl-vector-set-all *beta-hat* (/ *K*))

(define *Sigma-zero* (new-gsl-matrix *n* *n*)) ;;starting variances estimates
(gsl-matrix-set-identity *Sigma-zero*)
(define *Sigma-hat* (new-matrixArray *K*))
(matrixArray-setitem *Sigma-hat* 0 *Sigma-zero*)
(matrixArray-setitem *Sigma-hat* 1 *Sigma-zero*)

(define (make-target d S)
  "Build a gaussian-mixture target distribution
   d: distance between means
   S: ratio between variances"
  (let
      ((w (new-gsl-vector 2))
       (Sigma1 (new-gsl-matrix *n* *n*))
       (Sigma2 (new-gsl-matrix *n* *n*))
       (mu1 (new-gsl-vector *n*))
       (mu2 (new-gsl-vector *n*))
       (pi1 '())
       (pi2 '())
       (pi-mix '())
       (pi-array (new-mvnormArray 2)))
    (gsl-vector-set-all w 0.5)
    (gsl-matrix-set-identity Sigma1)
    (gsl-matrix-memcpy Sigma2 Sigma1)
    (gsl-matrix-scale Sigma2 S)
    (gsl-vector-set-all mu1 (- d))
    (gsl-vector-set-all mu2 d)
    (set! pi1 (new-mcmclib-mvnorm-lpdf mu1 (gsl-matrix-data-get Sigma1)))
    (set! pi2 (new-mcmclib-mvnorm-lpdf mu2 (gsl-matrix-data-get Sigma2)))
    (mvnormArray-setitem pi-array 0 pi1)
    (mvnormArray-setitem pi-array 1 pi2)
    (set! pi-mix (new-mcmclib-mixnorm-lpdf w pi-array))
    (lambda (x)
      (mcmclib-mixnorm-lpdf-compute pi-mix x))))

(define (simulate-amh d S sampler-builder)
  (let
      ((x (new-gsl-vector *n*)) ;;current chain state
       (target (make-target d S))
       (sampler '())
       (monitor '()))
    (gsl-vector-set-all x 0.0)
    (set! monitor (new-mcmclib-monitor x))
    (set! sampler (sampler-builder (mcmclib-guile-lpdf-cb) (guile-to-voidptr target)
                                   x))
    (do-ec (: i *T0*)
           (mcmclib-amh-update sampler))
    (do-ec (: i *N*)
           (begin
             (mcmclib-amh-update sampler)
             (mcmclib-monitor-update monitor)))
    monitor))

(define (simulate-raptor d S)
  (let*
      ((mu-hat (new-vectorArray *K*))
       (mu-vec (vector (new-gsl-vector *n*) (new-gsl-vector *n*)))
       (raptor-builder (lambda (fun fun-data x)
                         (mcmclib-raptor-alloc
                          *rng*
                          fun fun-data
                          x *T0*
                          *Sigma-zero*
                          *beta-hat*
                          mu-hat
                          *Sigma-hat*))))
    (gsl-vector-set-all (vector-ref mu-vec 0) (+ (* d -2) -1))
    (gsl-vector-set-all (vector-ref mu-vec 1) (+ (* d 2) 1))
    (do-ec (: i *K*)
           (vectorArray-setitem mu-hat i (vector-ref mu-vec i)))
    (simulate-amh d S raptor-builder)))

(define (simulate-gauss-am d S)
  (let*
      ((gauss-am-builder (lambda (fun fun-data x)
                        (mcmclib-gauss-am-alloc
                         *rng*
                         fun fun-data
                         x *Sigma-zero* *T0*))))
  (simulate-amh d S gauss-am-builder)))

(define (print-diags d S)
  "prints means, vars, ARs, MSJDs for RAPTOR and AM"
  (display "Gaussian AM") (newline)
  (mcmclib-monitor-fprintf-all (simulate-gauss-am d S) (stdout))
  (display "RAPTOR") (newline)
  (mcmclib-monitor-fprintf-all (simulate-raptor d S) (stdout)))

(print-diags 0 4)
(print-diags 3 1)
(print-diags 3 4)

(define (gv2v gv)
  "from gsl-vector to vector"
  (vector-ec (: i (gsl-vector-size-get gv)) (gsl-vector-get gv i)))

;;MSE comparison
(define (estimate-MSE replicas simulator)
  "estimate MSE of the sample mean estimator"
  (let
      ((mse (new-gsl-vector *n*))
       (mean-i (new-gsl-vector *n*))
       (update-sum-sq (lambda (curr-sum new-x)
                        (let
                            ((tmp-x (new-gsl-vector (gsl-vector-size-get new-x))))
                          (gsl-vector-memcpy tmp-x new-x)
                          (gsl-vector-mul tmp-x new-x)
                          (gsl-vector-add curr-sum tmp-x)))))
    (gsl-vector-set-all mse 0.0)
    (do-ec (: i replicas)
           (begin
             (mcmclib-monitor-get-means (simulator) mean-i)
             (update-sum-sq mse mean-i)))
    (gsl-vector-scale mse (/ replicas))
    (gv2v mse)))

(estimate-MSE 100 (lambda () (simulate-raptor 0 4)))
;;#(0.0468550464263089 0.051909728421788 0.0492880597329444 0.0376887694923336 0.0510109076993617)
(estimate-MSE 100 (lambda () (simulate-gauss-am 0 4)))
;;#(0.0417742622344229 0.0403734540469312 0.0437711141991529 0.0397311348579993 0.0435686305425603)

(estimate-MSE 100 (lambda () (simulate-raptor 3 4)))
(estimate-MSE 100 (lambda () (simulate-gauss-am 3 4)))

(estimate-MSE 100 (lambda () (simulate-raptor 3 1)))
(estimate-MSE 100 (lambda () (simulate-gauss-am 3 1)))
