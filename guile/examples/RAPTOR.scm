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
(define mu1 (v2gv #(0.2 0.2)))
(vectorArray-setitem *mu-hat* 0 mu1)
(define mu2 (v2gv #(0.8 0.8)))
(vectorArray-setitem *mu-hat* 1 mu2)

(define (ll2M ll)
  (let* ((n (length ll))
         (p (length (car ll)))
         (M (new-gsl-matrix n p)))
    (do-ec (: i n) (: j p)
           (gsl-matrix-set M i j (list-ref (list-ref ll i) j)))
    M))
(define *Sigma-zero* (ll2M '((0.01 0) (0 0.01))))
(matrixArray-setitem *Sigma-hat* 0 *Sigma-zero*)
(matrixArray-setitem *Sigma-hat* 1 *Sigma-zero*)

(define (dunif x)
  (let
      ((x0 (gsl-vector-get x 0))
       (x1 (gsl-vector-get x 0)))
    (if
     (and
      (and
       (>= x0 0)
       (<= x0 1))
      (and
       (>= x0 0)
       (<= x0 1)))
     0.0
     (log 0.0))))

(define *sampler* (mcmclib-raptor-alloc *rng*
                                        (mcmclib-guile-lpdf-cb) (guile-to-voidptr dunif)
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

(time (lambda ()
        (do-ec (: i 10000)
               (mcmclib-amh-update *sampler*))))
