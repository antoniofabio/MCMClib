;;
;;Compute the empirical CDF 'live' on a running chain
;;
(set! %load-path (cons "." %load-path))

(use-modules (srfi srfi-1)
	     (srfi srfi-42)
             (swig gsl)
             (swig mcmclib)
             (ice-9 format))

(define (gv2v gv)
  (vector-ec (: i (gsl-vector-size-get gv)) (gsl-vector-get gv i)))

(define (ecdf-init set)
  "init the ECDF data"
  (let*
      ((d (gsl-vector-size-get (vector-ref set 0))) ;domain dimension
       (Fn (new-gsl-vector (vector-length set))))
    (gsl-vector-set-all Fn 0.0)
    (vector set Fn 0)))

(define (ecdf-compute ecdf)
  "compute ECDF"
  (let*
      ((Fn (vector-ref ecdf 1))
       (ans (new-gsl-vector (gsl-vector-size-get Fn))))
    (gsl-vector-memcpy ans Fn)
    (gsl-vector-scale ans (/ (vector-ref ecdf 2)))
    (gv2v ans)))

(define (ecdf-update ecdf x-new)
  "update ecdf using point 'x-new'"
  (let*
      ((set (vector-ref ecdf 0))
       (Fn (vector-ref ecdf 1))
       (n (vector-ref ecdf 2))
       (zi (new-gsl-vector (gsl-vector-size-get (vector-ref set 0)))))
    (do-ec (: i (vector-length set))
           (begin
             (gsl-vector-memcpy zi (vector-ref set i))
             (gsl-vector-sub zi x-new)
             (gsl-vector-set Fn i (+ (gsl-vector-get Fn i) (if (= (gsl-vector-ispos zi) 0) 0 1)))))
    (vector set Fn (+ n 1))))

;;test it
(define (v2gv v)
  (let*
      ((n (vector-length v))
       (gv (new-gsl-vector n)))
    (do-ec (: i n)
           (gsl-vector-set gv i (vector-ref v i)))
    gv))

(define r (new-gsl-rng (gsl-rng-default)))
(define dim 3)
(define set (vector (v2gv #(0 0 0)) (v2gv #(0 1 1)) (v2gv #(0.2 1 1))
                    (v2gv #(0.8 1 1)) (v2gv #(1 1 1))))

(define (empirical-CDF X0 N sampler)
  (let
      ((ecdf (ecdf-init X0)))
    (do-ec (: i N)
           (set! ecdf (ecdf-update ecdf (sampler))))
    (ecdf-compute ecdf)))

(define ecdf  (empirical-CDF set 1000
               (lambda () (v2gv (vector-ec (: i dim) (gsl-rng-uniform r))))))
(define (pair i) (vector (vector-ref ecdf i) (gv2v (vector-ref set i))))
(do-ec (: i (vector-length set)) (format #t "~a\n" (pair i)))
