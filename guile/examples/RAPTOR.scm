(set! %load-path (cons "." %load-path))

(use-modules (srfi srfi-1)
	     (srfi srfi-42)
             (swig gsl)
             (swig mcmclib))

(define *n* 2) ;;target distribution dimension
(define *K* 2) ;;number of mixture components
(define *rng* (new-gsl-rng (gsl-rng-default)))
(define *x* (new-gsl-vector *n*)) ;;current chain state
(define *T0* 10000) ;;burn-in length

(define *beta-hat* (new-gsl-vector *K*))
(define *mu-hat* (new-vectorArray *K*))
(define *Sigma-hat* (new-matrixArray *K*))

(gsl-vector-set-all *beta-hat* (/ *K*))
(define (v2gv v)
  (let*
      ((n (vector-length v))
       (gv (new-gsl-vector n)))
    (do-ec (: i n)
           (gsl-vector-set gv i (vector-ref v i)))
    gv))
(define mu1 (v2gv #(-5 -5)))
(vectorArray-setitem *mu-hat* 0 mu1)
(define mu2 (v2gv #(5 5)))
(vectorArray-setitem *mu-hat* 1 mu2)

(define (ll2M ll)
  (let* ((n (length ll))
         (p (length (car ll)))
         (M (new-gsl-matrix n p)))
    (do-ec (: i n) (: j p)
           (gsl-matrix-set M i j (list-ref (list-ref ll i) j)))
    M))
(define *Sigma-zero* (ll2M '((1 0) (0 1))))
(matrixArray-setitem *Sigma-hat* 0 *Sigma-zero*)
(matrixArray-setitem *Sigma-hat* 1 *Sigma-zero*)

;;Build a gaussian-mixture target distribution
;;d: distance between means
;;S: ratio between variances
(define (make-target d S)
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

(define target (make-target 0 4))

(gsl-vector-set-all *x* 0.5)

;;Make a RAPTOR sampler for distance 'd', ratio 'S'
(define *sampler* (mcmclib-raptor-alloc *rng*
                                        (mcmclib-guile-lpdf-cb) (guile-to-voidptr target)
                                        *x* *T0*
                                        *Sigma-zero*
                                        *beta-hat*
                                        *mu-hat*
                                        *Sigma-hat*))

(define (time fun)
  (let*
      ((start (current-time))
       (thrash (fun))
       (end (current-time)))
    (- end start)))

(do-ec (: i *T0*)
       (mcmclib-amh-update *sampler*))

(define mon (new-mcmclib-monitor *x*))
(time (lambda ()
        (do-ec (: i 100000)
               (begin
                 (mcmclib-amh-update *sampler*)
                 (mcmclib-monitor-update mon)))))
(mcmclib-monitor-fprintf-all mon (stdout)) ;;prints means, vars, ARs, MSJDs
