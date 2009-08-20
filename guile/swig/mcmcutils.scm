(define-module (swig mcmcutils))

(use-modules (srfi srfi-1)
	     (srfi srfi-42)
             (swig gsl)
             (swig mcmclib)
             (oop goops))

(export v2gv gv2v gM2M M2gM va2ca ma2ca)

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
(define (M2gM M)
  "convert a scheme matrix into a gsl matrix"
  (let*
      ((nrows (vector-length M))
       (ncols (vector-length (vector-ref M 0)))
       (gM (new-gsl-matrix nrows ncols)))
    (do-ec (: i nrows)
           (let
               ((row (vector-ref M i)))
             (do-ec (: j ncols)
                    (gsl-matrix-set gM i j (vector-ref row j)))))
    gM))

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

(define (xor a b)
  (or (and a (not b))
      (and (not a) b)))
(define (square x) (* x x))
(define (dist x y)
  "euclidean distance between two gsl vectors"
  (let*
      ((n (gsl-vector-size-get x))
       (z (new-gsl-vector n)))
    (gsl-vector-memcpy z x)
    (gsl-vector-sub z y)
    (do-ec (: i n) (gsl-vector-set z i (square (gsl-vector-get z i))))
    (sum-ec (: i n) (gsl-vector-get z i))))
(define (get-dists x pts)
  "compute distances between point x and list of points 'pts'. gsl vectors"
  (map (lambda (y) (dist x y)) pts))
(define (which-min x)
  "0-based index of the minimum value in list 'x'"
  (let
      ((m (apply min x)))
    (list-index (lambda (y) (= y m)) x)))

;;
;;keep references of referenced objects to avoid premature garbage collection
;;
(define-class <swig-obj> () (c-ref #:init-keyword #:c-ref #:getter get-c-ref))
(define-class <distrfun> ()
  (fun-ptr #:init-keyword #:fun-ptr)
  (fun-data-ptr #:init-keyword #:fun-data-ptr)
  (keep #:init-keyword #:keep))
(define (make-guile-distrfun fun)
  (make <distrfun>
    #:fun-ptr (mcmclib-guile-lpdf-cb)
    #:fun-data-ptr (guile-to-voidptr fun)
    #:keep fun))
(define-class <mh> (<swig-obj>)
  (rng #:init-keyword #:rng)
  (distrfun #:init-keyword #:distrfun)
  (x #:init-keyword #:x))
(define-method (update (obj <mh>)) (mcmclib-mh-update (get-c-ref obj)))
(define-class <amh> (<mh>))
(define-method (update (obj <amh>)) (mcmclib-amh-update (get-c-ref obj)))
(export <swig-obj> <distrfun> <mh> <amh> make-guile-distrfun update get-c-ref)

(use-syntax (ice-9 syncase))

(define (symbol-concatenate lst)
  (string->symbol (string-concatenate (map symbol->string lst))))
(export symbol-concatenate)

(define-syntax make-mh
  (syntax-rules ()
    ((make-mh sub-type rng-in distrfun-obj-in x-in rest ...)
     (let
         ((constructor-name (symbol-concatenate (list 'mcmclib- sub-type '-alloc)))
          (rng rng-in)
          (x x-in)
          (distrfun-obj distrfun-obj-in))
     (make <mh>
       #:c-ref ((eval constructor-name (interaction-environment))
                rng
                (slot-ref distrfun-obj 'fun-ptr)
                (slot-ref distrfun-obj 'fun-data-ptr)
                x rest ...)
       #:rng rng
       #:distrfun distrfun-obj
       #:x x)))))
(export-syntax make-mh)

(define-syntax make-amh
  (syntax-rules ()
    ((make-amh sub-type rng-in distrfun-obj-in x-in rest ...)
     (let
         ((constructor-name (symbol-concatenate (list 'mcmclib- sub-type '-alloc)))
          (rng rng-in)
          (x x-in)
          (distrfun-obj distrfun-obj-in))
     (make <amh>
       #:c-ref ((eval constructor-name (interaction-environment))
                rng
                (slot-ref distrfun-obj 'fun-ptr)
                (slot-ref distrfun-obj 'fun-data-ptr)
                x rest ...)
       #:rng rng
       #:distrfun distrfun-obj
       #:x x)))))
;;use as follows:
;;
;;(make-amh 'raptor rng (make-guile-distrfun (lambda (x) 0.0) x
;;          t0 Sigma_zero beta_hat mu_hat Sigma_hat)
;;(make-amh 'gauss-am rng distrfun distrfun-data x sigma_zero t0)
(export-syntax make-amh)

;;
;;wrap monitor objects
;;
(define-class <monitor> (<swig-obj>) (x #:init-keyword #:x))
(define-method (update (obj <monitor>))
  (mcmclib-monitor-update (get-c-ref obj)))
(define-class <monitor-ecdf> (<swig-obj>) (x #:init-keyword #:x))
(define-method (update (obj <monitor-ecdf>))
  (mcmclib-monitor-ecdf-update (get-c-ref obj) (slot-ref obj 'x)))
(export <monitor> <monitor-ecdf>)

;;
;;transition frequencies of a finite state machine
;;
(define-class <transition-matrix> () counts state)
(define (make-transition-matrix n init)
  "builds a new transition matrix. 'n' is the number of states, 'init' the initial value"
  (let
      ((ans (make <transition-matrix>)))
    (slot-set! ans 'state init)
    (slot-set! ans 'counts (make-array 0 n n))
    ans))
(define-method (update (obj <transition-matrix>) new-state)
  (let*
      ((old-state (slot-ref obj 'state))
       (old-count (array-ref (slot-ref obj 'counts) old-state new-state)))
    (array-set! (slot-ref obj 'counts) (+ old-count 1) old-state new-state)
    (slot-set! obj 'state new-state)
    obj))
