(define (make-mvnorm mu Sigma)
  "build a multivariate normal distrib. fun."
  (let
      ((c-data (new-mcmclib-mvnorm-lpdf mu (gsl-matrix-data-get Sigma))))
    (lambda (x)
      (mcmclib-mvnorm-lpdf-compute c-data x))))

(define (make-mixnorm beta mu Sigma)
  "build a normal mixture distrib."
  (let
      ((beta-list (vector->list beta))
       (funs (list-ec
              (: i (vector-length beta))
              (make-mvnorm (vector-ref mu i) (vector-ref Sigma i)))))
    (lambda (x)
      (log (reduce + 0
                   (map
                    (lambda (fun coef)
                      (*
                       (exp (fun x))
                       coef))
                    funs
                    beta-list))))))

;;Monitoring functions
(define (which-region x pts-list)
  (which-min (get-dists x pts-list)))
(define-class <monitor-region> () pts-list x trans-mat)
(define (make-monitor-region pts-list x)
  (let*
      ((n (length pts-list))
       (trans-mat (make-transition-matrix n (which-region x pts-list)))
       (ans (make <monitor-region>)))
    (slot-set! ans 'pts-list pts-list)
    (slot-set! ans 'trans-mat trans-mat)
    (slot-set! ans 'x x)
    ans))
(define-method (update (obj <monitor-region>))
  (update (slot-ref obj 'trans-mat) (which-region (slot-ref obj 'x) (slot-ref obj 'pts-list))))

(define (update-gv-mean old-mean new-value old-n)
  (gsl-vector-scale old-mean old-n)
  (gsl-vector-add old-mean new-value)
  (gsl-vector-scale old-mean (/ (+ old-n 1.0))))
(define (update-gM-mean old-mean new-value old-n)
  (gsl-matrix-scale old-mean old-n)
  (gsl-matrix-add old-mean new-value)
  (gsl-matrix-scale old-mean (/ (+ old-n 1.0))))

(define-class <chain-diags> () vectors matrices n)
(define (make-chain-diags dim matrices)
  (let ((ans (make <chain-diags>)))
    (slot-set! ans 'vectors (list-ec (: i 3) (make-filled-vector 0.0 dim)))
    (slot-set! ans 'matrices matrices)
    (slot-set! ans 'n 0.0)
    ans))
(define-method (update (obj <chain-diags>) new-diags)
  (let
      ((n (slot-ref obj 'n)))
    (map (lambda (x y) (update-gv-mean x (v2gv y) n))
         (slot-ref obj 'vectors)
         (car new-diags))
    (map (lambda (x y) (update-gM-mean x (M2gM y) n))
         (slot-ref obj 'matrices)
         (cadr new-diags))
    (slot-set! obj 'n (+ n 1.0))
    obj))
(define-method (write (obj <chain-diags>) port)
  (write
   (list
    (map gv2v (slot-ref obj 'vectors))
    (map gM2M (slot-ref obj 'matrices)))
   port))
