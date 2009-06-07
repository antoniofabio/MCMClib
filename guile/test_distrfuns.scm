(load "main.scm")

(define rng (new-gsl-rng (gsl-rng-default)))
(define v (new-gsl-vector 2))
(gsl-vector-set-all v 0.0)
(define S (new-gsl-matrix 2 2))
(gsl-matrix-set-identity S)
(define p (new-mcmclib-mvnorm-lpdf v (gsl-matrix-data-get S)))
(mcmclib-mvnorm-lpdf-compute p v)

(mcmclib-mvnorm rng S v)

;;MCAR distribution
(define n 95)
(define p 3)
(define W (new-gsl-matrix n n))
(gsl-matrix-set-all W 0.0)
(do-ec (: i (- n 1))
       (begin
         (gsl-matrix-set W i (+ i 1) 1.0)
         (gsl-matrix-set W (+ i 1) i 1.0)))
(define lpdf-data (new-mcmclib-mcar-tilde-lpdf p W))
(define e (new-gsl-vector (* p n)))
(gsl-vector-set-all e 0.0)
(mcmclib-mcar-tilde-lpdf-compute lpdf-data e)

(define mod (new-mcmclib-mcar-model lpdf-data e))
(define as (mcmclib-mcar-tilde-lpdf-alphasigmag-get lpdf-data))
(define a12s (mcmclib-mcar-tilde-lpdf-alpha12sigma-get lpdf-data))
(define y (new-gsl-vector (* p p)))
(gsl-vector-set-all y 0.0)
(do-ec (: i p)
       (gsl-vector-set y
                       (- (* p p) (+ i 1))
                       (+ -2 (+ (* (/ i p) 0.9) 0.05))))
(gsl-vector-memcpy a12s y)
(gsl-vector-set-all as 0.0)
(do-ec (: i p)
       (gsl-vector-set as
                       (-
                        (+
                         (* p
                            (/ (- p 1)
                               2))
                         p)
                        (+ i 1))
                       (+ i 0.05)))

(mcmclib-mcar-model-alpha12sigma-lpdf mod y)
