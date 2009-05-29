;(use-modules (ice-9 debugging traps)
;	     (ice-9 gds-client)
;	     (ice-9 debugging example-fns)
;	     (srfi srfi-42)
;             (swig mcmclib))

;(do-ec (: i 3) (display i))

;(install-trap (make <procedure-trap>
;		#:behaviour gds-debug-trap
;		#:procedure fact1))

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
(define *phi-par-samplers* (make-vector 2))
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

(define *util-lambdas* (make-vector *n*)) ;this is needed to prevent garbage collection
(define (make-phij-fcond j)
  (vector-set! *util-lambdas* j
               (lambda (phij)
                 (mcmclib-mcar-model-phi-fcond *mcar-model* j phij)))
  (makeVoidPtr (vector-ref *util-lambdas* j)))

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

(define *X* (new-gsl-matrix *np* *p*))
(gsl-matrix-set-zero *X*)
(do-ec (: i *n*) (: j *p*)
       (gsl-matrix-set *X*
                       (+ j (* i *p*))
                       j
                       1.0))

(define *y* (new-gsl-vector *np*))
(gsl-vector-fscanf (fopen "y_2.dat" "r") *y*)
(define *offset* (new-gsl-vector *np*))
(gsl-vector-fscanf (fopen "offset_2.dat" "r") *offset*)
(do-ec (: i *np*)
       (gsl-vector-set *offset* i (log (gsl-vector-get *offset* i))))

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
  (do-ec (: i *n*)
         (let ((v (new-gsl-vector *p*)))
           (gsl-vector-set-all v -1.0)
           (vector-set! *phi-i-samplers* i
                        (am-sampler (guile-distrfun)
                                    (make-phij-fcond i)
                                    v))
           (vector-set! *phi-i-vecs* i v)))
  (set! *beta-sampler*
        (mcmclib-pmodel-sampler-sampler-get
         (new-mcmclib-pmodel-sampler *X* *y* *offset* *rng* 1e-3 *T0*))))

(define (update-phi)
  (do-ec (: i *n*)
         (let ((smp (vector-ref *phi-i-samplers* i))
               (v (vector-ref *phi-i-vecs* i)))
           (mcmclib-amh-update smp)
           (gsl-copy-subvec *phi* v (* i *p*))
           (update-offset))))

(define (update-phi-pars)
  (do-ec (: i 2)
         (mcmclib-amh-update (vector-ref *phi-par-samplers* i))))

(define (update-beta) (mcmclib-amh-update *beta-sampler*))

(define tmp (list-ec (: i 10)
                     (begin
                       (update-phi-pars)
                       (gsl-vector-get *phi* 0))))

(define (update-all N)
  (do-ec (: i N)
         (begin
           (update-phi-pars)
           ;;    (update-phi)
           (update-beta))))
(define (main)
  (update-all *N*))

(init-chains)
(define st (current-time))
(main)
(define en (current-time))
(display "elapsed time (seconds):")(newline)
(display (- en st))(newline)
