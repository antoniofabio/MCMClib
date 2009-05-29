(set! %load-path (cons "." %load-path))

(use-modules (srfi srfi-1)
	     (srfi srfi-42)
             (swig gsl)
             (swig mcmclib))

;;list -> g-vector
(define (l2v ll)
  (let*
      ((n (length ll))
       (v (new-gsl-vector n)))
    (do-ec (: i n) (gsl-vector-set v i (list-ref ll i)))
    v))

;;g-matrix -> g-vector
(define (M2v M)
  (let* ((n (gsl-matrix-size1-get M))
         (p (gsl-matrix-size2-get M))
         (v (new-gsl-vector (* n p))))
    (do-ec (: i n) (: j p)
           (gsl-vector-set v
                        (+ (* i p) j)
                        (gsl-matrix-get M i j)))
    v))

;;list of lists -> matrix
(define (ll2M ll)
  (let* ((n (length ll))
         (p (length (car ll)))
         (M (new-gsl-matrix n p)))
    (do-ec (: i n) (: j p)
           (gsl-matrix-set M i j (list-ref (list-ref ll i) j)))
    M))
(ll2M '((1 2 3) (4 5 6)))

;;g-vector -> list
(define (v2l v)
  (list-ec (: i (gsl-vector-size-get v))
           (gsl-vector-get v i)))
(v2l (l2v '(3 1 4)))

;;get matrix row as a list
(define (M2l M i)
  (list-ec (: j (gsl-matrix-size2-get M))
           (gsl-matrix-get M i j)))
(define M (ll2M '((1 2 3) (4 5 6))))
(M2l M 0)
(M2l M 1)

;;g-matrix -> list of lists
(define (M2ll M)
  (list-ec (: i (gsl-matrix-size1-get M))
           (list-ec (: j (gsl-matrix-size2-get M))
                    (gsl-matrix-get M i j))))
(M2ll (ll2M '((1 2 3) (4 5 6))))

;;display g-vector
(define (dv v) (display (v2ll v)) (newline))
;;display g-Matrix
(define (dM M) (display (M2ll M)) (newline))

(use-modules (swig mcmclib))

(define P 3)
(define Psi (new-gsl-matrix P P))
(gsl-matrix-set-identity Psi)
(define p (mcmclib-iwishart-lpdf-alloc Psi P))

(define v (new-gsl-vector (* P P)))
(gsl-vector-set-all v 1.0)
(mcmclib-iwishart-lpdf-compute p v)

(define (gsl-copy-subvec dest src offset)
  (do-ec (: i (gsl-vector-size-get src))
         (gsl-vector-set dest (+ offset i) (gsl-vector-get src i))))
