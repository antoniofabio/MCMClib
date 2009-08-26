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

;;Update mean of a list of objects on which and '(add x y)'
;;  and a '(scale x s)' method is defined.
(define (update-means old-values new-values old-n)
  (map (lambda (x y) (update-mean x y old-n)) old-values new-values))
