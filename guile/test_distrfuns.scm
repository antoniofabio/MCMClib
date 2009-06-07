(load "main.scm")

(define rng (new-gsl-rng (gsl-rng-default)))
(define v (new-gsl-vector 2))
(gsl-vector-set-all v 0.0)
(define S (new-gsl-matrix 2 2))
(gsl-matrix-set-identity S)
(define p (new-mcmclib-mvnorm-lpdf v (gsl-matrix-data-get S)))
(mcmclib-mvnorm-lpdf-compute p v)

(mcmclib-mvnorm rng S v)
