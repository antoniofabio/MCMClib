(define-module (swig mcmcutils))

(use-modules (srfi srfi-1)
	     (srfi srfi-42)
             (swig gsl-utils)
             (swig mcmclib)
             (oop goops)
	     (ice-9 optargs)
             (ice-9 pretty-print))
(use-modules ((swig gsl) :renamer gsl-renamer))

(use-syntax (ice-9 syncase))

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
(define-method (mul (x <vector>) . rest)
  (let-optional rest
   ((y x))
   (vector-ec (: i (vector-length x)) (* (vector-ref x i) (vector-ref y i)))))
(export add scale mul)

(define-public (update-mean old-value new-data old-n)
  "update mean value 'old-value' based on sample size 'old-n'
   using the new data point 'new-data'"
  (scale (add (scale old-value old-n) new-data) (/ (+ 1.0 old-n))))

(define-public (xor a b)
  (or (and a (not b))
      (and (not a) b)))

(define-public (square x) (* x x))

(define-public (get-dists x pts)
  "compute distances between point x and list of points 'pts'. gsl vectors"
  (map (lambda (y) (gv-dist x y)) pts))

(define-public (which-min x)
  "0-based index of the minimum value in list 'x'"
  (let
      ((m (apply min x)))
    (list-index (lambda (y) (= y m)) x)))

;;
;;keep references of referenced objects to avoid premature garbage collection
;;
(define-class <swig-obj> ()
  (c-ref #:init-keyword #:c-ref #:getter get-c-ref)
  (subtype #:init-keyword #:subtype #:getter get-subtype))
(define-method (free (obj <swig-obj>))
  (let ((destructor-name (symbol-concatenate (list (get-subtype obj) '-free))))
    ((eval destructor-name (interaction-environment)) (get-c-ref obj))))

(define-public (symbol-concatenate lst)
  (string->symbol (string-concatenate (map symbol->string lst))))

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
       (d (gv-size-get (mcmclib-monitor-x-get cobj)))
       (v (new-gv d))
       (fun-name (symbol-concatenate (list 'mcmclib-monitor-get- what))))
    ((eval fun-name (interaction-environment)) cobj v)
    (gv2v v)))
(define-method (write (obj <monitor>) port)
  (let*
      ((cobj (get-c-ref obj))
       (d (gv-size-get (mcmclib-monitor-x-get cobj)))
       (v (new-gv d)))
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
