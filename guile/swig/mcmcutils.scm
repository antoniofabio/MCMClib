(define-module (swig mcmcutils))

(use-modules (srfi srfi-1)
	     (srfi srfi-42)
             (swig gsl)
             (swig mcmclib)
             (oop goops))

(export v2gv gv2v gM2M va2ca ma2ca)

(define (v2gv v)
  "convert the scheme vector 'v' into a gsl vector"
  (let*
      ((n (vector-length v))
       (gv (new-gsl-vector n)))
    (do-ec (: i n)
           (gsl-vector-set gv i (vector-ref v i)))
    gv))
(define (gv2v gv)
  "convert the gsl vector 'gv' into a scheme vector"
  (vector-ec (: i (gsl-vector-size-get gv)) (gsl-vector-get gv i)))

(define (gM2M gM)
  "convert a gsl matrix into a scheme matrix"
  (vector-ec (: i (gsl-matrix-size1-get gM))
             (vector-ec (: j (gsl-matrix-size2-get gM))
                        (gsl-matrix-get gM i j))))

(define (va2ca va)
  "convert a vector of g-vectors into a C array of g-vectors"
  (let*
      ((n (vector-length va))
       (ca (new-vectorArray n)))
    (do-ec (: i n)
           (vectorArray-setitem ca i (vector-ref va i)))
    ca))
(define (ma2ca ma)
  "convert a vector of g-matrices into a C array of g-matrices"
  (let*
      ((n (vector-length ma))
       (ca (new-matrixArray n)))
    (do-ec (: i n)
           (matrixArray-setitem ca i (vector-ref ma i)))
    ca))

;;
;;keep references of referenced objects to avoid premature garbage collection
;;
(define-class <swig-obj> () (c-ref #:init-keyword #:c-ref #:getter get-c-ref))
(define-class <amh> (<swig-obj>) (rng #:init-keyword #:rng) (x #:init-keyword #:x))
(define-method (update (obj <amh>)) (mcmclib-amh-update (get-c-ref obj)))
(export <swig-obj> <amh> update)

(use-syntax (ice-9 syncase))

(define (symbol-concatenate lst)
  (string->symbol (string-concatenate (map symbol->string lst))))
(export symbol-concatenate)

(define-syntax make-amh
  (syntax-rules ()
    ((make-amh sub-type rng distrfun distrfun-data x ...)
     (let
         ((constructor-name (symbol-concatenate (list 'mcmclib- sub-type '-alloc))))
     (make <amh>
       #:c-obj ((eval constructor-name (interaction-environment))
                rng distrfun distrfun-data x ...)
       #:rng rng
       #:x (car (list x ...)))))))
;;use as follows:
;;
;;(make-amh 'raptor rng distrfun distrfun-data x
;;          t0 Sigma_zero beta_hat mu_hat Sigma_hat)
;;(make-amh 'gauss-am rng distrfun distrfun-data x sigma_zero t0)
(export-syntax make-amh)
