(load-extension "./libschememcmclib.so" "SWIG_init")

(define P 3)
(define Psi (gsl-matrix-alloc P P))
(gsl-matrix-set-identity Psi)
(define m 3)
(define p (mcmclib-iwishart-lpdf-alloc Psi m))

;;list -> g-vector
(define (ll2v ll)
  (let
      ((n (length ll))
       (v (gsl-vector-alloc n))
    (let loop ((i 0))
      (gsl-vector-set v i (list-ref ll i))
      (if (i < (- n 1)) (loop (+ i 1)))))))

(define (M2v M)
  (let ((n (gsl-matrix-size1-get M))
        (v (gsl-vector-alloc (* n n))))
    (let loop-i ((i 0))
      (let loop-j ((j 0))
        (gsl-vector-set v
                        (+ (* i n) j)
                        (gsl-matrix-get M i j))
        (if (< j (- n 1)) (loop-j (+ j 1))))
      (if (< i (- n 1)) (loop-i (+ i 1))))))

(define v (gsl-vector-alloc (* P P)))
(gsl-vector-set-all v 1.0)
(mcmclib-iwishart-lpdf-compute p v)
