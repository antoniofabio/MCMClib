(load "main.scm")

(define DIM 2)
(define N 5)
(define K 2)

(define .mu (vector-ec (: i K) (v2gv #( 1 1))))
(define .sigma (vector-ec (: i K)
                          (let
                              ((M (new-gsl-matrix 2 2)))
                            (gsl-matrix-set-zero M)
                            M)))
(define mu-hat (new-vectorArray K))
(define Sigma-hat (new-matrixArray K))
(do-ec (: i K)
       (begin
         (vectorArray-setitem mu-hat i (vector-ref .mu i))
         (matrixArray-setitem Sigma-hat i (vector-ref .sigma i))))
(define w-hat (new-gsl-vector K))
(gsl-vector-set-all w-hat (/ K))
(define m (new-mcmclib-mixem-online mu-hat Sigma-hat w-hat 0.5 2))
(define y (v2gv #(1 1)))
;this halts with a numerical error:
;(do-ec (: i N)
;       (mcmclib-mixem-online-update m y))
