;;
;;Compute the empirical CDF 'live' on a running chain
;;
(set! %load-path (cons "." %load-path))

(use-modules (srfi srfi-1)
	     (srfi srfi-42)
             (swig gsl)
             (swig mcmclib))

(define (gv2v gv)
  (vector-ec (: i (gsl-vector-size-get gv)) (gsl-vector-get gv i)))

(define (ecdf-init set)
  "init the ECDF data"
  (let*
      ((d (gsl-vector-size-get (vector-ref set 0))) ;domain dimension
       (Fn (new-gsl-vector (vector-length set))))
    (gsl-vector-set-all Fn 0.0)
    (vector set Fn 0)))

(define (ecdf-compute ecdf x)
  "compute ECDF on point 'x'"
  (gv2v (gsl-vector-scale (vector-ref ecdf 1) (vector-ref ecdf 2))))

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
             (gsl-vector-add (gsl-vector-get Fn i) (if (gsl-vector-ispos zi) 0 1))))
    (vector set Fn (+ n 1))))
