;;
;;compute mixture of normals ECDF
;;
(set! %load-path (cons "." %load-path))

(use-modules (srfi srfi-1)
	     (srfi srfi-42)
             (swig gsl)
             (swig mcmclib)
             (ice-9 format))

(define *n* 5) ;;target distribution dimension
(define *K* 2) ;;number of mixture components
(define *T0* 100) ;;burn-in length
(define *N* 1000) ;;number of iterations after burn-in
(define *rng* (new-gsl-rng (gsl-rng-default))) ;;random number generator

(load "empiricalCDF.scm")

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

(define (make-sampler d S sampler-builder)
  (let*
      ((target (make-target d S))
       (x (new-gsl-vector *n*))
       (sampler-obj (sampler-builder (mcmclib-guile-lpdf-cb) (guile-to-voidptr target) x)))
    (lambda ()
      (mcmclib-amh-update sampler-obj)
      x)))

(define (make-raptor-sampler d S)
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
    (make-sampler d S raptor-builder)))

(define (make-gauss-am-sampler d S)
  (let*
      ((gauss-am-builder (lambda (fun fun-data x)
                        (mcmclib-gauss-am-alloc
                         *rng*
                         fun fun-data
                         x *Sigma-zero* *T0*))))
  (make-sampler d S gauss-am-builder)))

(define B 1000)
(define (generate-sample)
  (v2gv (vector-ec (: i *n*) (gsl-rng-uniform *rng*))))
(define X0 (vector-ec (: i B) (generate-sample)))
(define ecdf-true (empirical-CDF X0 1000 (make-gauss-am-sampler 1 1)))
(define ecdf-RAPTOR (empirical-CDF X0 100 (make-raptor-sampler 1 1)))
(define ecdf-gauss-AM (empirical-CDF X0 100 (make-gauss-am-sampler 1 1)))
(define (compute-dist a b)
  (let*
      ((ga (v2gv a))
       (gb (v2gv b))
       (n (vector-length a))
       (ans 0.0))
    (do-ec (: i n)
           (begin
             (set! ans
                   (+ ans
                      (abs
                       (-
                        (gsl-vector-get ga i)
                        (gsl-vector-get gb i)))))))
    (/ ans n)))
(define (get-distance ecdf-estimator-fun)
  (compute-dist ecdf-true (ecdf-estimator-fun)))
(define raptor-distances
  (vector-ec (: i 10)
             (get-distance (lambda ()
                             (empirical-CDF X0 100 (make-raptor-sampler 1 1))))))
(define gauss-am-distances
  (vector-ec (: i 10)
             (get-distance (lambda ()
                             (empirical-CDF X0 100 (make-gauss-am-sampler 1 1))))))
