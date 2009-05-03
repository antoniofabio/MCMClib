(load-extension "./libschememcmclib.so" "SWIG_init")

;;list -> g-vector
(define (l2v ll)
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

;;list of lists -> matrix
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

(define (dec n) (- n 1))

;;generate a sequence
(define (seq i)
  (cond
   ((= i 1) (list 1))
   (else
    (append
     (seq (- i 1))
     (list i)))))

;;g-vector -> list
(define (v2l v)
  (map
   (lambda (i) (gsl-vector-get v i))
   (map dec (seq (gsl-vector-size-get v)))))
(v2l (l2v '(3 1 4)))

;;get matrix row as a list
(define (M2l M i)
  (map
   (lambda (j) (gsl-matrix-get M i j))
   (map dec (seq (gsl-matrix-size2-get M)))))
(define M (ll2M '((1 2 3) (4 5 6))))
(M2l M 0)
(M2l M 1)

;;g-matrix -> list of lists
(define (M2ll M)
  (map
   (lambda (i) (M2l M i))
   (map dec (seq (gsl-matrix-size1-get M)))))
(M2ll (ll2M '((1 2 3) (4 5 6))))

;;display g-vector
(define (dv v) (display (v2ll v)) (newline))
;;display g-Matrix
(define (dM M) (display (M2ll M)) (newline))

(define P 3)
(define Psi (gsl-matrix-alloc P P))
(gsl-matrix-set-identity Psi)
(define m P)
(define p (mcmclib-iwishart-lpdf-alloc Psi m))

(define v (gsl-vector-alloc (* P P)))
(gsl-vector-set-all v 1.0)
(mcmclib-iwishart-lpdf-compute p v)
