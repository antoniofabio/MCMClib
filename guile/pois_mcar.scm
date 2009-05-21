(load "main.scm")

;;
;;Check Poisson modeling
;;
(define *N* 1000)
(define *THIN* 10)
(define *T0* 1000)
(define *V0* 0.4)
(define *p* 3)
(define *n* 95)
(define *np* (* *n* *p*))

(define *rng* (new-gsl-rng (gsl-rng-default)))
(define *phi-par-samplers* (list))
(define *phi-i-samplers* (make-vector *n*))
(define *phi-i-vecs* (make-vector *n*))
(define *beta-sampler* (list))
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
  (makeVoidPtr
   (lambda (phij)
     (mcmclib-mcar-model-phi-fcond *mcar-model* j phij))))

(define *alpha12sigma* (mcmclib-mcar-tilde-lpdf-alpha12sigma-get *mcar-lik*))
(define *alphasigmag* (mcmclib-mcar-tilde-lpdf-alphasigmag-get *mcar-lik*))

(define (Sigma0 dim)
  (let ((ans (new-gsl-matrix dim dim)))
    (gsl-matrix-set-identity ans)
    (gsl-matrix-scale ans (/ *V0* dim))
    ans))

(define (am-sampler lpdf lpdf-data x)
  (mcmclib-gauss-am-alloc *rng* lpdf lpdf-data x
                          (Sigma0 (gsl-vector-size-get x))
                          *T0*))

(define *X* (new-gsl-matrix (* *n* *p*) *p*))
(gsl-matrix-set-zero *X*)
(do ((i 0 (+ i *p*))) ((>= i (* *n* *p*)))
  (do ((j 0 (1+ j))) ((>= j *p*))
    (gsl-matrix-set *X* (+ i j) j 1.0)))

(define *y* (new-gsl-vector (* *n* *p*)))
(gsl-vector-fscanf (fopen "y_2.dat" "r") *y*)
(define *offset* (new-gsl-vector (* *n* *p*)))
(gsl-vector-fscanf (fopen "offset_2.dat" "r") *offset*)
(do ((i 0 (1+ i))) ((>= i (* *n* *p*)))
  (gsl-vector-set *offset* i (log (gsl-vector-get *offset* i))))

(define (gsl-copy-subvec dest src offset)
  (do ((i 0 (1+ i))) ((>= i (gsl-vector-size-get src)))
    (gsl-vector-set dest (+ offset i) (gsl-vector-get src i))))

(define (init-chains)
  (gsl-vector-set-all *alpha12sigma* -1.0)
  (gsl-vector-set-all *alphasigmag* -1.0)
  (set! *phi-par-samplers*
        (list
         (am-sampler (mcmclib-mcar-model-alpha12sigma-lpdf-cb)
                     *mcar-model*
                     *alpha12sigma*)
         (am-sampler (mcmclib-mcar-model-alphasigma-lpdf-cb)
                     *mcar-model*
                     *alphasigmag*)))
  (do ((i 0 (1+ i))) ((>= i *n*))
    (vector-set! *phi-i-vecs* i (new-gsl-vector *p*))
    (vector-set! *phi-i-samplers* i
                 (am-sampler (guile-distrfun)
                             (make-phij-fcond i)
                             (vector-ref *phi-i-vecs* i))))
  (set! *beta-sampler*
        (mcmclib-pmodel-sampler-sampler-get
         (mcmclib-pmodel-sampler-alloc *X* *y* *offset* *rng* 1e-3 *T0*))))

(define (update-phi)
  (do ((i 0 (1+ i))) ((>= i))
    (mcmclib-amh-update (vector-ref *phi-i-samplers* i))
    (gsl-copy-subvec *phi* (vector-ref *phi-i-vecs* i))
    (update-offset)))

(define (update-phi-pars)
  (map
   (lambda (smp) (mcmclib-amh-update smp) (mcmclib-mcar-tilde-lpdf-update-vcov *mcar-lik*))
   *phi-par-samplers*))

(define (update-beta) (mcmclib-amh-update *beta-sampler*))

(define (main)
  (init-chains)
  (do ((i 0 (1+ i))) ((>= i *N*))
    (update-phi-pars)
    (update-phi)
    (update-beta)))

(define st (current-time))
(main)
(define en (current-time))
(display "elapsed time (seconds):")(newline)
(display (- en st))
