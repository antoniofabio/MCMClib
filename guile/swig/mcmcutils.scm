(define-module (swig mcmcutils))

(use-modules (srfi srfi-1)
	     (srfi srfi-42)
             (swig gsl)
             (swig mcmclib))

(export v2gv gv2v gM2M va2ca ma2ca)

(define (v2gv v)
  "convert the scheme vector 'v' into a gsl vector"
  (let*
      ((n (vector-length v))
       (gv (new-gsl-vector n)))
    (do-ec (: i n)
           (gsl-vector-set gv i (vector-ref v i)))
    gv))
(define (gv2v gv)
  "convert the gsl vector 'gv' into a scheme vector"
  (vector-ec (: i (gsl-vector-size-get gv)) (gsl-vector-get gv i)))

(define (gM2M gM)
  "convert a gsl matrix into a scheme matrix"
  (vector-ec (: i (gsl-matrix-size1-get gM))
             (vector-ec (: j (gsl-matrix-size2-get gM))
                        (gsl-matrix-get gM i j))))

(define (va2ca va)
  "convert a vector of g-vectors into a C array of g-vectors"
  (let*
      ((n (vector-length va))
       (ca (new-vectorArray n)))
    (do-ec (: i n)
           (vectorArray-setitem ca i (vector-ref va i)))
    ca))
(define (ma2ca ma)
  "convert a vector of g-matrices into a C array of g-matrices"
  (let*
      ((n (vector-length ma))
       (ca (new-matrixArray n)))
    (do-ec (: i n)
           (matrixArray-setitem ca i (vector-ref ma i)))
    ca))
