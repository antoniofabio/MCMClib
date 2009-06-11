(load "main.scm")

;;
;;Check Poisson modeling
;;
(define *T0* 1000)
(define *V0* 0.4)
(define *p* 3)
(define *n* 95)
(define *np* (* *n* *p*))

;;Load data
(define (with-file fname mode fun)
  (let*
      ((f (fopen fname mode))
       (ans (fun f)))
    (fclose f)
    ans))
(define (file->vector fname gvec)
  (with-file fname "r"
             (lambda (f)
               (gsl-vector-fscanf f gvec))))

(define *offset* (new-gsl-vector *np*))
(define *denom* (new-gsl-vector *np*))
(file->vector "offset_2.dat" *denom*)
(define *y* (new-gsl-vector *np*))
(file->vector "y_2.dat" *y*)

(define *W* (new-gsl-matrix *n* *n*))
(with-file "W_2.dat" "r"
           (lambda (f)
             (gsl-matrix-fscanf f *W*)))

;;Setup design matrix
(define *X* (new-gsl-matrix *np* *p*))
(gsl-matrix-set-zero *X*)
(do-ec (: i *n*) (: j *p*)
       (gsl-matrix-set *X*
                       (+ j (* i *p*))
                       j
                       1.0))

(do-ec (: i *np*)
       (gsl-vector-set *denom* i (log (gsl-vector-get *denom* i))))

;;Helper functions
(define (update-offset)
  (gsl-vector-memcpy *offset* *denom*)
  (gsl-vector-add *offset* *phi*))

;;Init RNG
(define *rng* (new-gsl-rng (gsl-rng-default)))

;;Setup Gamma and B distribs
(define *phi* (new-gsl-vector *np*))
(define *mcar-lik* (new-mcmclib-mcar-tilde-lpdf *p* *W*))
(define *mcar-model* (new-mcmclib-mcar-model *mcar-lik* *phi*))

(define *alpha12sigma* (mcmclib-mcar-tilde-lpdf-alpha12sigma-get *mcar-lik*))
(define *alphasigmag* (mcmclib-mcar-tilde-lpdf-alphasigmag-get *mcar-lik*))
(define *phi-par-samplers* (list))

;;Setup phi-i full-conditionals
(define (make-fcond i)
  (lambda (phi-i)
    (mcmclib-mcar-model-phi-fcond *mcar-model* i phi-i)))
(define *fconds*
  (vector-ec (: i *n*) (make-fcond i)))

;;phi-i samplers builder
(define Si (new-gsl-matrix *p* *p*))
(gsl-matrix-set-identity Si)
(gsl-matrix-scale Si *V0*)
(define *phi-samplers* (list))
(define (make-phi-sampler i)
  (let*
      ((phi-i (new-gsl-vector *p*))
       (thrash (gsl-vector-set-all phi-i -1))
       (mh (mcmclib-gauss-mrw-alloc
            *rng*
            (mcmclib-guile-lpdf-cb)
            (guile-to-voidptr (vector-ref *fconds* i))
            phi-i
            Si))
       (ip (* i *p*)))
    (lambda ()
      (mcmclib-mh-update mh)
      (gsl-copy-subvec *phi* phi-i ip)
      (update-offset))))

;;Setup beta sampler
(define *beta-save* (list))
(define *beta-sampler* (list))

;;Generic AM sampler builder
(define (Sigma0 dim)
  (let ((ans (new-gsl-matrix dim dim)))
    (gsl-matrix-set-identity ans)
    (gsl-matrix-scale ans (/ *V0* dim))
    ans))
(define (am-sampler lpdf lpdf-data x)
  (mcmclib-gauss-am-alloc *rng* lpdf lpdf-data x
                          (Sigma0 (gsl-vector-size-get x))
                          *T0*))

;;Init chains: starting values and samplers
(define (init-chains)
  (gsl-vector-set-all *alpha12sigma* -1.0)
  (gsl-vector-set-all *alphasigmag* -1.0)
  (set! *phi-par-samplers*
        (vector
         (am-sampler (mcmclib-mcar-model-alpha12sigma-lpdf-cb)
                     *mcar-model*
                     *alpha12sigma*)
         (am-sampler (mcmclib-mcar-model-alphasigma-lpdf-cb)
                     *mcar-model*
                     *alphasigmag*)))
  (set! *phi-samplers*
        (vector-ec (: i *n*)
                   (make-phi-sampler i)))
;  (set! *beta-save* (new-mcmclib-pmodel-sampler *X* *y* *offset* *rng* 1e-3 *T0*))
;  (set! *beta-sampler*
;        (mcmclib-pmodel-sampler-sampler-get *beta-save*))
  )

;;Updating functions
(define (update-phi)
  (do-ec (: i *n*)
         ((vector-ref *phi-samplers* i))))
(define (update-phi-pars)
  (do-ec (: i 2)
         (mcmclib-amh-update (vector-ref *phi-par-samplers* i))))
(define (update-beta) (mcmclib-amh-update *beta-sampler*))

(define (update-all N)
  (do-ec (: i N)
         (begin
           (update-phi-pars)
           (update-phi)
;           (update-beta)
           )))

(define (see-phi)
  (vector-ec (: i 10)
             (gsl-vector-get *phi* i)))

(init-chains)
(update-all 1000)
