(load "main.scm")

;;
;;Check Poisson modeling
;;
(define *N* 50000)
(define *THIN* 10)
(define *T0* 1000)
(define *V0* 0.4)
(define *p* 3)
(define *n* 95)

(define *samplers* (list))

(define *W* (new-gsl-matrix *n* *n*))
(define *mcar-lik* (new-mcmclib-mcar-tilde-lpdf *p* *W*))
(define *phi* (new-gsl-vector (* *p* *n*)))
(define *mcar-model* (new-mcmclib-mcar-model *mcar-lik* *phi*))
