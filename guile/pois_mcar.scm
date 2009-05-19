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
(define *np* (* *n* *p*))

(define *samplers* (list))
(define *W* (new-gsl-matrix *n* *n*))
(define *mcar-lik* (new-mcmclib-mcar-tilde-lpdf *p* *W*))
(define *phi* (new-gsl-vector *np*))
(define *mcar-model* (new-mcmclib-mcar-model *mcar-lik* *phi*))

(define *denom* (new-gsl-vector *np*))
(define *offset* (new-gsl-vector *np*))

(define (update-offset)
  (gsl-vector-memcpy *offset* *denom*)
  (gsl-vector-add *offset* *phi*))

(define (make-phij-fcond j)
  (lambda (phij)
    (mcmclib_mcar_model_phi_fcond *mcar_model* j phij)))

(define (init-chains))
(define (free-chains))

(define (main))
