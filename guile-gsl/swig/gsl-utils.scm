(define-module (swig gsl-utils)
  :use-syntax (ice-9 optargs))
(use-modules (srfi srfi-42)
	     (ice-9 syncase)
	     (ice-9 optargs)
	     (swig gsl))

(define-public (sv->gv sv)
  "convert a scheme vector into a gsl vector"
  (let*
      ((size (vector-length sv))
       (ans (new-gsl-vector size)))
    (do-ec (:range i size) (gsl-vector-set ans i (vector-ref sv i)))
    ans))
(define-public (gv->sv gv)
  "convert a gsl vector into a scheme vector"
  (vector-ec (:range i (gsl-vector-size-get gv)) (gsl-vector-get gv i)))

;; gsl vector eager comprehension
(define-syntax gv-ec
  (syntax-rules (nested)
    ((gv-ec expression)
     (gv-ec (nested) expression))

    ((gv-ec qualifier expression)
     (sv->gv (vector-ec qualifier expression)))))
(export-syntax gv-ec)

;; gsl vector generator
(define-syntax :gv
  (syntax-rules (index)
    ((:gv cc var (index var1) arg)
     (:vector cc var (index var1) (gv->sv arg)))
    ((:gv cc var arg)
     (:vector cc var (gv->sv arg)))))
(export-syntax :gv)

;; gsl vector index generator
(define-syntax :gv-along
  (syntax-rules ()
    ((:gv-along cc var arg)
     (:range cc var (gsl-vector-size-get arg)))))
(export-syntax :gv-along)

(define-public (gv-fold x0 x f)
  (fold-ec x0 (:gv xi x) xi f))
(define-public (gv-clone gv) (gv-ec (:gv x gv) x))
(define-public (gv-sum gv) (gv-fold gv 0 +))
(define-public (gv-min gv) (min-ec (:gv x gv) x))
(define-public (gv-max gv) (max-ec (:gv x gv) x))
(define-public (gv-map f gv) (gv-ec (:gv x gv) (f x)))
(define-public (gv-map-ip f gv)
  "map the function 'f' over the gsl vector 'gv' in place"
  (do-ec (:gv gvi (index i) gv)
	 (gsl-vector-set gv i (f gvi)))
  gv)

(define-public (gv-copy-subvec dest src offset)
  "copy a gvector into a subset of another gvector"
  (do-ec (: i (gsl-vector-size-get src))
         (gsl-vector-set dest (+ offset i) (gsl-vector-get src i))))

(define-public (gv-dist x y)
  "euclidean distance between two gsl vectors"
  (sum-ec (:parallel (:gv xi x) (:gv yi y)) (:let di (- xi yi))
	  (* di di)))

(define-public (gm->sm gm)
  "convert a gsl matrix into a scheme matrix"
  (let*
      ((size1 (gsl-matrix-size1-get gm))
       (size2 (gsl-matrix-size2-get gm))
       (ans (make-array 0.0 size1 size2)))
    (do-ec (:range i size1) (:range j size2)
	   (array-set! ans (gsl-matrix-get gm i j) i j))
    ans))
(define-public (sm->gm sm)
  "convert a scheme matrix into a gsl matrix"
  (let*
      ((sizes (array-dimensions sm))
       (size1 (car sizes))
       (size2 (cadr sizes))
       (ans (new-gsl-matrix 0.0 size1 size2)))
    (do-ec (:range i size1) (:range j size2)
	   (gsl-matrix-set ans i j (array-ref sm i j)))
    ans))

(define-public (gm-row gm i)
  "clone row 'i' of gsl matrix 'gm'"
  (gv-ec (:range j (gsl-matrix-size2-get gm)) (gsl-matrix-get gm i j)))
(define-public (gm-rows gm)
  "clone rows of matrix 'gm' into a scheme vector of gsl vectors"
  (vector-ec (:range i (gsl-matrix-size1-get gm)) (gm-row gm i)))
;; gsl matrix rows generator
(define-syntax :gm-rows
  (syntax-rules (index)
    ((:gm-rows cc var (index var1) arg)
     (:vector cc var (index var1) (gm-rows arg)))
    ((:gm-rows cc var arg)
     (:vector cc var (gm-rows arg)))))
(export-syntax :gm-rows)

(define-public (gm-col gm j)
  "clone column 'i' of gsl matrix 'gm'"
  (gv-ec (:range i (gsl-matrix-size1-get gm)) (gsl-matrix-get gm i j)))
(define-public (gm-cols gm)
  "clone columns of matrix 'gm' into a scheme vector of gsl vectors"
  (vector-ec (:range j (gsl-matrix-size2-get gm)) (gm-col gm j)))
;; gsl matrix columns generator
(define-syntax :gm-cols
  (syntax-rules (index)
    ((:gm-cols cc var (index var1) arg)
     (:vector cc var (index var1) (gm-rows arg)))
    ((:gm-cols cc var arg)
     (:vector cc var (gm-rows arg)))))
(export-syntax :gm-cols)

(define-public (gm-memcpy-gv gm gv)
  "copy contents of gsl vector 'gv' into gsl matrix 'gm'"
  (let
      ((size1 (gsl-matrix-size1-get gm))
       (size2 (gsl-matrix-size2-get gm)))
    (do-ec (:range i size1) (:range j size2)
	   (gsl-matrix-set gm i j (gsl-vector-get gv (+ (* j size1) i))))
    gm))

(define-public (gm-clone gm)
  (let*
      ((size1 (gsl-matrix-size1-get gm))
       (size2 (gsl-matrix-size2-get gm))
       (ans (new-gsl-matrix size1 size2)))
    (do-ec (:range i size1) (:range j size2)
	   (gsl-matrix-set ans i j (gsl-matrix-get gm i j)))
    ans))
(define-public (gm-map-ip gm f)
  (do-ec (:range i (gsl-matrix-size1-get gm))
	 (:range j (gsl-matrix-size2-get gm))
	 (gsl-matrix-set gm i j (f (gsl-matrix-get gm i j))))
  gm)
(define-public (gm-map gm f)
  (gm-map-ip (gm-clone gm) f))

(define-public (make-gm-diag dim . rest)
  "make a diagonal gsl matrix of dimension 'dim' and optional value 'value'"
  (let-optional rest
   ((value 1.0))
   (let
       ((ans (new-gsl-matrix dim dim)))
     (gsl-matrix-set-all ans 0.0)
     (do-ec (:range i dim) (gsl-matrix-set ans i i value))
     ans)))
