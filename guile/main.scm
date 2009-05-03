(load-extension "./libschememcmclib.so" "SWIG_init")

(define (ll2v ll)
  (let*
      ((n (length ll))
       (v (gsl-vector-alloc n)))
       (display n) (display v) (newline)))

;;list -> g-vector
(define (ll2v ll)
  (let*
      ((n (length ll))
       (v (gsl-vector-alloc n)))
    (let loop ((i 0))
      (gsl-vector-set v i (list-ref ll i))
      (if (< i (- n 1)) (loop (+ i 1))))
    v))

;;g-matrix -> g-vector
(define (M2v M)
  (let* ((n (gsl-matrix-size1-get M))
         (p (gsl-matrix-size2-get M))
         (v (gsl-vector-alloc (* n p))))
    (let loop-i ((i 0))
      (let loop-j ((j 0))
        (gsl-vector-set v
                        (+ (* i p) j)
                        (gsl-matrix-get M i j))
        (if (< j (- n 1)) (loop-j (+ j 1))))
      (if (< i (- p 1)) (loop-i (+ i 1))))
    v))

;;list -> matrix
(define (ll2M ll)
  (let* ((n (length ll))
         (p (length (car ll)))
         (M (gsl-matrix-alloc n p)))
    (let loop-i ((i 0))
      (let loop-j ((j 0))
        (let ((row (list-ref ll i)))
          (gsl-matrix-set M i j (list-ref row j)))
        (if (< j (- p 1)) (loop-j (+ j 1))))
      (if (< i (- n 1)) (loop-i (+ i 1))))
    M))
(ll2M '((1 2 3) (4 5 6)))

;;display g-vector
(define (dv v) (display (v2ll v)))
;;display g-Matrix
(define (dM M) (display (M2ll M)))

(define P 3)
(define Psi (gsl-matrix-alloc P P))
(gsl-matrix-set-identity Psi)
(define m P)
(define p (mcmclib-iwishart-lpdf-alloc Psi m))

(define v (gsl-vector-alloc (* P P)))
(gsl-vector-set-all v 1.0)
(mcmclib-iwishart-lpdf-compute p v)
