(set! %load-path (cons "." %load-path))

(use-modules (srfi srfi-1)
	     (srfi srfi-42)
             (swig gsl)
             (swig mcmclib))

(define *n* 2) ;;target distribution dimension
(define *K* 2) ;;number of mixture components
(define *T0* 10000) ;;burn-in length
(define *N* 100000) ;;number of iterations after burn-in

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

(define (simulate-raptor d S)
  (let
      (
       (rng (new-gsl-rng (gsl-rng-default))) ;;random number generator
       (x (new-gsl-vector *n*)) ;;current chain state
       (mu-hat (new-vectorArray *K*))
       (mu-vec (vector (new-gsl-vector *n*) (new-gsl-vector *n*)))
       (target (make-target d S))
       (sampler '())
       (monitor '())
       )
    (gsl-vector-set-all x 0.0)
    (set! monitor (new-mcmclib-monitor x))
    (gsl-vector-set-all (vector-ref mu-vec 0) (+ (* d -2) -1))
    (gsl-vector-set-all (vector-ref mu-vec 1) (+ (* d 2) 1))
    (do-ec (: i *K*)
           (vectorArray-setitem mu-hat i (vector-ref mu-vec i)))
    (set! sampler (mcmclib-raptor-alloc rng
                                        (mcmclib-guile-lpdf-cb) (guile-to-voidptr target)
                                        x *T0*
                                        *Sigma-zero*
                                        *beta-hat*
                                        mu-hat
                                        *Sigma-hat*))
    (do-ec (: i *T0*)
           (mcmclib-amh-update sampler))
    (do-ec (: i *N*)
           (begin
             (mcmclib-amh-update sampler)
             (mcmclib-monitor-update monitor)))
    monitor
  ))

(define (time fun)
  (let*
      ((start (current-time))
       (thrash (fun))
       (end (current-time)))
    (- end start)))

(define *monitor* (simulate-raptor 0 4))
(mcmclib-monitor-fprintf-all *monitor* (stdout)) ;;prints means, vars, ARs, MSJDs
