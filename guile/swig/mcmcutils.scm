(define-module (swig mcmcutils))

(use-modules (srfi srfi-1)
	     (srfi srfi-42)
             (swig gsl)
             (swig mcmclib)
             (oop goops)
             (ice-9 pretty-print))
(use-syntax (ice-9 syncase))

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
  (let*
      ((nr (gsl-matrix-size1-get gM))
       (nc (gsl-matrix-size2-get gM))
       (ans (make-array 0.0 nr nc)))
    (do-ec (: i nr) (: j nc)
           (array-set! ans (gsl-matrix-get gM i j) i j))
    ans))
(define-public (matrix-numrows M)
  "get the number of rows of matrix M"
  (car (array-dimensions M)))
(define-public (matrix-numcolumns M)
  "get the number of columns of matrix M"
  (cadr (array-dimensions M)))
(define (M2gM M)
  "convert a scheme matrix into a gsl matrix"
  (let*
      ((nr (matrix-numrows M))
       (nc (matrix-numcolumns M))
       (gM (new-gsl-matrix nr nc)))
    (do-ec (: i nr) (: j nc)
           (gsl-matrix-set gM i j (array-ref M i j)))
    gM))

(define-public (matrix-add a b)
  (let
      ((ans (M2gM a)))
    (gsl-matrix-add ans (M2gM b))
    (gM2M ans)))
(define-public (matrix-scale a b)
  (let
      ((ans (M2gM a)))
    (gsl-matrix-scale ans b)
    (gM2M ans)))

(define-method (add (x <vector>) (y <vector>))
  (vector-ec (: i (vector-length x)) (+ (vector-ref x i) (vector-ref y i))))
(define-method (scale (x <vector>) s)
  (vector-ec (: i (vector-length x)) (* (vector-ref x i) s)))
(define-method (add (x <array>) (y <array>)) (matrix-add x y))
(define-method (scale (x <array>) s) (matrix-scale x s))
(define-method (add (a <list>) (b <list>)) (map add a b))
(define-method (scale (a <list>) s) (map (lambda (x) (scale x s)) a))
(define-method (add (a <number>) (b <number>)) (+ a b))
(define-method (scale (a <number>) s) (* a s))
(export add scale)

(define-public (update-mean old-value new-data old-n)
  "update mean value 'old-value' based on sample size 'old-n'
   using the new data point 'new-data'"
  (scale (add (scale old-value old-n) new-data) (/ (+ 1.0 old-n))))

(define (gsl-copy-subvec dest src offset)
  "copy a gvector into a subset of another gvector"
  (do-ec (: i (gsl-vector-size-get src))
         (gsl-vector-set dest (+ offset i) (gsl-vector-get src i))))

(define (va2ca va)
  "convert a vector of g-vectors into a C array of g-vectors"
  (let*
      ((n (vector-length va))
       (ca (new-vectorArray n)))
    (do-ec (: i n)
           (vectorArray-setitem ca i (vector-ref va i)))
    ca))
(define (ca2va ca size)
  (vector-ec (: i size) (gv2v (vectorArray-getitem ca i))))
(define (ma2ca ma)
  "convert a vector of g-matrices into a C array of g-matrices"
  (let*
      ((n (vector-length ma))
       (ca (new-matrixArray n)))
    (do-ec (: i n)
           (matrixArray-setitem ca i (vector-ref ma i)))
    ca))
(define (ca2ma ca size)
  (vector-ec (: i size) (gM2M (matrixArray-getitem ca i))))

(define (xor a b)
  (or (and a (not b))
      (and (not a) b)))

(define (square x) (* x x))

(define (dist x y)
  "euclidean distance between two gsl vectors"
  (sum-ec (: i (gsl-vector-size-get x))
          (square (- (gsl-vector-get x i) (gsl-vector-get y i)))))

(define (get-dists x pts)
  "compute distances between point x and list of points 'pts'. gsl vectors"
  (map (lambda (y) (dist x y)) pts))

(define (which-min x)
  "0-based index of the minimum value in list 'x'"
  (let
      ((m (apply min x)))
    (list-index (lambda (y) (= y m)) x)))

(define (diag value size)
  "gsl diagonal matrix 'size x size'"
  (let
      ((ans (new-gsl-matrix size size)))
    (gsl-matrix-set-identity ans)
    (gsl-matrix-scale ans value)
    ans))

(define (make-filled-vector value dim)
  "make gsl vector of dimension 'dim' filled with 'value'"
  (let
      ((ans (new-gsl-vector dim)))
    (gsl-vector-set-all ans value)
    ans))

(export xor square dist get-dists which-min diag make-filled-vector)

;;
;;keep references of referenced objects to avoid premature garbage collection
;;
(define-class <swig-obj> ()
  (c-ref #:init-keyword #:c-ref #:getter get-c-ref)
  (subtype #:init-keyword #:subtype #:getter get-subtype))
(define-method (free (obj <swig-obj>))
  (let ((destructor-name (symbol-concatenate (list (get-subtype obj) '-free))))
    ((eval destructor-name (interaction-environment)) (get-c-ref obj))))
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
(export <swig-obj> <distrfun> <mh> <amh> make-guile-distrfun update
        get-c-ref get-subtype free)

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
       #:subtype (symbol-concatenate (list 'mcmclib- sub-type))
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
       #:subtype (symbol-concatenate (list 'mcmclib- sub-type))
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
(define (make-monitor x)
  (make <monitor> #:x x #:c-ref (new-mcmclib-monitor x)))
(define-method (update (obj <monitor>))
  (mcmclib-monitor-update (get-c-ref obj)))
(define-method (get-monitor-value (obj <monitor>) what)
  "get values 'what' (a symbol) from monitor object 'obj' as a scheme vector"
  (let*
      ((cobj (get-c-ref obj))
       (d (gsl-vector-size-get (mcmclib-monitor-x-get cobj)))
       (v (new-gsl-vector d))
       (fun-name (symbol-concatenate (list 'mcmclib-monitor-get- what))))
    ((eval fun-name (interaction-environment)) cobj v)
    (gv2v v)))
(define-method (write (obj <monitor>) port)
  (let*
      ((cobj (get-c-ref obj))
       (d (gsl-vector-size-get (mcmclib-monitor-x-get cobj)))
       (v (new-gsl-vector d)))
    (format #t "<monitor>\n")
    (format #t "monitoring: ~a\n" (slot-ref obj 'x))
    (format #t "means:\n")
    (pretty-print (get-monitor-value obj 'means))
    (format #t "variances:\n")
    (pretty-print (get-monitor-value obj 'vars))
    (format #t "acceptance rates:\n")
    (pretty-print (get-monitor-value obj 'ar))
    (format #t "mean squared jumping distances:\n")
    (pretty-print (get-monitor-value obj 'msjd))))

(define-class <monitor-ecdf> (<swig-obj>) (x #:init-keyword #:x))
(define-method (update (obj <monitor-ecdf>))
  (mcmclib-monitor-ecdf-update (get-c-ref obj) (slot-ref obj 'x)))
(export <monitor> get-monitor-value <monitor-ecdf> make-monitor)

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
(export <transition-matrix> make-transition-matrix)
