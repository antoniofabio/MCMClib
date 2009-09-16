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
(define (v2gv v)
  (let*
      ((n (vector-length v))
       (gv (new-gsl-vector n)))
    (do-ec (: i n)
           (gsl-vector-set gv i (vector-ref v i)))
    gv))


(define *N-X0* 1000)
(define *dim* 5)
(define *X0* (new-gsl-matrix *N-X0* *dim*))
(let
    ((f (fopen "mixnorm_iid_sample.dat" "r")))
  (gsl-matrix-fscanf f *X0*)
  (fclose f))
;;compute CDF
(define *monitor* (new-mcmclib-monitor-ecdf *X0*))
(do-ec (: i *N-X0*)
       (let*
           ((xi (v2gv (vector-ec (: j *dim*) (gsl-matrix-get *X0* i j)))))
         (mcmclib-monitor-ecdf-update *monitor* xi)))
(define *Fn-tgt* (mcmclib-monitor-ecdf-Fn-get *monitor*))
