(load "main.scm")

;;
;;Example MCAR distribution, with GRW sampling
;;
(define *N* 1000)
(define *THIN* 10)
(define *V0* 0.05)
(define *p* 6)
(define *n* 5)
(define *np* (* *n* *p*))

(define *rng* (new-gsl-rng (gsl-rng-default)))
(define *W* (new-gsl-matrix *n* *n*))
(define *mcar-lik* (new-mcmclib-mcar-tilde-lpdf *p* *W*))
(define *phi* (new-gsl-vector *np*))
(gsl-matrix-fscanf (fopen "mcar_grw_W.dat" "r") *W*)
(gsl-vector-fscanf (fopen "mcar_grw_phi.dat" "r") *phi*)
(define *mcar-model* (new-mcmclib-mcar-model *mcar-lik* *phi*))
(define *samplers* (make-vector 2))

(define *alpha12sigma* (mcmclib-mcar-tilde-lpdf-alpha12sigma-get *mcar-lik*))
(define *alphasigmag* (mcmclib-mcar-tilde-lpdf-alphasigmag-get *mcar-lik*))
(define *monitors*
  (vector
   (new-mcmclib-monitor *alpha12sigma*)
   (new-mcmclib-monitor *alphasigmag*)))

(define (Sigma0 dim)
  (let ((ans (new-gsl-matrix dim dim)))
    (gsl-matrix-set-identity ans)
    (gsl-matrix-scale ans (/ *V0* dim))
    ans))

(define (sampler lpdf lpdf-data x)
  (mcmclib-gauss-mrw-alloc *rng* lpdf lpdf-data x
                           (Sigma0 (gsl-vector-size-get x))))

(define (init-chains)
  (gsl-vector-set-all *alpha12sigma* -1.0)
  (gsl-vector-set-all *alphasigmag* -1.0)
  (set! *samplers*
        (vector
         (sampler (mcmclib-mcar-model-alpha12sigma-lpdf-cb)
                     *mcar-model*
                     *alpha12sigma*)
         (sampler (mcmclib-mcar-model-alphasigma-lpdf-cb)
                  *mcar-model*
                  *alphasigmag*))))

(define (update-all N)
  (do-ec (: n N)
         (do-ec (: i 2)
                (begin
                  (mcmclib-mh-update (vector-ref *samplers* i))
                  (mcmclib-monitor-update (vector-ref *monitors* i))
                  ))))

(define (main)
  (update-all *N*))

(init-chains)
(define st (current-time))
(main)
(define en (current-time))
(display "elapsed time (seconds):")(newline)
(display (- en st))(newline)
