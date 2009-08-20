(define (xor a b)
  (or (and a (not b))
      (and (not a) b)))
(define (square x) (* x x))
(define (dist x y)
  (let*
      ((n (gsl-vector-size-get x))
       (z (new-gsl-vector n)))
    (gsl-vector-memcpy z x)
    (gsl-vector-sub z y)
    (do-ec (: i n) (gsl-vector-set z i (square (gsl-vector-get z i))))
    (sum-ec (: i n) (gsl-vector-get z i))))
(define (get-dists x pts)
  (map (lambda (y) (dist x y)) pts))
(define (which-min x)
  (let
      ((m (apply min x)))
    (list-index (lambda (y) (= y m)) x)))

;transition frequencies of a finite state machine
(define-class <transition-matrix> () counts state)
(define (make-transition-matrix n init)
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
