(define-module (swig mcmcutils))

(use-modules (srfi srfi-1)
	     (srfi srfi-42)
             (swig gsl)
             (swig mcmclib)
             (oop goops)
             (ice-9 pretty-print))
(use-syntax (ice-9 syncase))

(define-public (v2gv v)
  "convert the scheme vector 'v' into a gsl vector"
  (let*
      ((n (vector-length v))
       (gv (new-gsl-vector n)))
    (do-ec (: i n)
           (gsl-vector-set gv i (vector-ref v i)))
    gv))
(define-public (gv2v gv)
  "convert the gsl vector 'gv' into a scheme vector"
  (vector-ec (: i (gsl-vector-size-get gv)) (gsl-vector-get gv i)))

(define-public (gM2M gM)
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
(define-public (M2gM M)
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

(define-public (gsl-copy-subvec dest src offset)
  "copy a gvector into a subset of another gvector"
  (do-ec (: i (gsl-vector-size-get src))
         (gsl-vector-set dest (+ offset i) (gsl-vector-get src i))))

(define-public (gM-clone x)
  (let
      ((y (new-gsl-matrix (gsl-matrix-size1-get x) (gsl-matrix-size2-get x))))
    (gsl-matrix-memcpy y x)
    y))

(define-public (gv-clone x)
  (let
      ((y (new-gsl-vector (gsl-vector-size-get x))))
    (gsl-vector-memcpy y x)
    y))

(define-public (gv-op-ip x f)
  "map function f to gsl-vector x, in place"
  (do-ec (: i (gsl-vector-size-get x))
         (gsl-vector-set x i (f (gsl-vector-get x i))))
  x)

(define-public (gv-op x f) (gv-op-ip (gv-clone x) f))

(define-public (gv-fold x x0 f)
  (fold-ec x0 (: i (gsl-vector-size-get x)) (gsl-vector-get x i) f))

(define-public (gv-sum x)
  "sum elements of vector 'x'"
  (gv-fold x 0 +))

(define-public (xor a b)
  (or (and a (not b))
      (and (not a) b)))

(define-public (square x) (* x x))

(define-public (dist x y)
  "euclidean distance between two gsl vectors"
  (sum-ec (: i (gsl-vector-size-get x))
          (square (- (gsl-vector-get x i) (gsl-vector-get y i)))))

(define-public (get-dists x pts)
  "compute distances between point x and list of points 'pts'. gsl vectors"
  (map (lambda (y) (dist x y)) pts))

(define-public (which-min x)
  "0-based index of the minimum value in list 'x'"
  (let
      ((m (apply min x)))
    (list-index (lambda (y) (= y m)) x)))

(define-public (diag value size)
  "gsl diagonal matrix 'size x size'"
  (let
      ((ans (new-gsl-matrix size size)))
    (gsl-matrix-set-identity ans)
    (gsl-matrix-scale ans value)
    ans))

(define-public (make-filled-vector value dim)
  "make gsl vector of dimension 'dim' filled with 'value'"
  (let
      ((ans (new-gsl-vector dim)))
    (gsl-vector-set-all ans value)
    ans))

;;
;;keep references of referenced objects to avoid premature garbage collection
;;
(define-class <swig-obj> ()
  (c-ref #:init-keyword #:c-ref #:getter get-c-ref)
  (subtype #:init-keyword #:subtype #:getter get-subtype))
(define-method (free (obj <swig-obj>))
  (let ((destructor-name (symbol-concatenate (list (get-subtype obj) '-free))))
    ((eval destructor-name (interaction-environment)) (get-c-ref obj))))
(define-class <mh> (<swig-obj>)
  (rng #:init-keyword #:rng)
  (distrfun #:init-keyword #:distrfun)
  (x #:init-keyword #:x))
(define-method (update (obj <mh>)) (mcmclib-mh-update (get-c-ref obj)))
(define-class <amh> (<mh>))
(define-method (update (obj <amh>)) (mcmclib-amh-update (get-c-ref obj)))
(export <swig-obj> <mh> <amh> update get-c-ref get-subtype free)

(define-public (symbol-concatenate lst)
  (string->symbol (string-concatenate (map symbol->string lst))))

(define-syntax make-mh
  (syntax-rules ()
    ((make-mh sub-type rng-in distrfun-in x-in rest ...)
     (let
         ((constructor-name (symbol-concatenate (list 'mcmclib- sub-type '-alloc)))
          (rng rng-in)
          (x x-in)
          (distrfun distrfun-in))
     (make <mh>
       #:c-ref ((eval constructor-name (interaction-environment))
                rng
                distrfun
                x rest ...)
       #:subtype (symbol-concatenate (list 'mcmclib- sub-type))
       #:rng rng
       #:distrfun distrfun
       #:x x)))))
(export-syntax make-mh)

(define-syntax make-amh
  (syntax-rules ()
    ((make-amh sub-type rng-in distrfun-in x-in rest ...)
     (let
         ((constructor-name (symbol-concatenate (list 'mcmclib- sub-type '-alloc)))
          (rng rng-in)
          (x x-in)
          (distrfun distrfun-in))
     (make <amh>
       #:c-ref ((eval constructor-name (interaction-environment))
                rng
                distrfun
                x rest ...)
       #:subtype (symbol-concatenate (list 'mcmclib- sub-type))
       #:rng rng
       #:distrfun distrfun
       #:x x)))))
;;use as follows:
;;
;;(make-amh 'gauss-am rng lambda (x) 0.0) x sigma_zero t0)
(export-syntax make-amh)

(define-public (make-rapt rng distrfun x t0 sigma-whole K sigma-local region-fun)
  "alloc a new 'rapt' sampler object"
  (make-amh 'rapt rng distrfun x t0 sigma-whole K sigma-local region-fun))

;;
;;wrap monitor objects
;;
(define-class <monitor> (<swig-obj>) (x #:init-keyword #:x))
(define-public (make-monitor x)
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
(export <monitor> get-monitor-value <monitor-ecdf>)

;;
;;transition frequencies of a finite state machine
;;
(define-class <transition-matrix> () counts state)
(define-public (make-transition-matrix n init)
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
(export <transition-matrix>)
