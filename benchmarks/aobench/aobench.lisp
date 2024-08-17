;; aobench.lisp
;;
;; Path tracing
;; Benchmark for numeric compuations, aggregate access
;;
;; Based on Gauche Scheme code, which is a port of original C code
;; by Syoyo Fujita and licensed under New BSD License.
;; https://github.com/syoyo/aobench

(cl:in-package #:benchmark-aobench)

(cl:in-package #:benchmark-aobench/native)

(cl:declaim (cl:optimize (cl:speed 3) (cl:safety 1)))

;;-------------------------------------------------
;; Common utilities
;;

(coalton-toplevel

(define-struct (Vec3 :t)
  (velts (arr:LispArray :t)))                   ;length 3

(define-struct (Mat3 :t)
  (melts (arr:LispArray :t)))                   ;length 9

;; convenience constructor
(declare v3 (types:RuntimeRepr :t => :t -> :t -> :t -> (Vec3 :t)))
(define (v3 x y z)
  (let ((type (types:runtime-repr (types:proxy-of x))))
    (Vec3 (lisp (arr:LispArray :t) (type x y z)
            (cl:make-array '(3) :initial-contents (cl:list x y z)
                                :element-type type)))))

(declare vset! ((Vec3 :t) -> :t -> :t -> :t -> Unit))
(define (vset! v x y z)
  (arr:set! (.velts v) 0 x)
  (arr:set! (.velts v) 1 y)
  (arr:set! (.velts v) 2 z))

(declare vcopy! ((Vec3 :t) -> (Vec3 :t) -> Unit))
(define (vcopy! v s)
  (arr:set! (.velts v) 0 (arr:aref (.velts s) 0))
  (arr:set! (.velts v) 1 (arr:aref (.velts s) 1))
  (arr:set! (.velts v) 2 (arr:aref (.velts s) 2)))

(declare vx ((Vec3 :t) -> :t))
(declare vy ((Vec3 :t) -> :t))
(declare vz ((Vec3 :t) -> :t))
(define (vx v) (arr:aref (.velts v) 0))
(define (vy v) (arr:aref (.velts v) 1))
(define (vz v) (arr:aref (.velts v) 2))

(declare vdup (types:RuntimeRepr :t => (Vec3 :t) -> (Vec3 :t)))
(define (vdup v)
  (v3 (vx v) (vy v) (vz v)))

(declare vdot (Num :t => (Vec3 :t) -> (Vec3 :t) -> :t))
(define (vdot v w)
  (+ (* (vx v) (vx w))
     (+ (* (vy v) (vy w))
        (* (vz v) (vz w)))))

(declare vscale ((types:RuntimeRepr :t) (Num :t) => (Vec3 :t) -> :t -> (Vec3 :t)))
(define (vscale v a)
  (v3 (* a (vx v))
      (* a (vy v))
      (* a (vz v))))

(declare vscale! ((Num :t) => (Vec3 :t) -> :t -> (Vec3 :t)))
(define (vscale! v a)
  (arr:set! (.velts v) 0 (* (arr:aref (.velts v) 0) a))
  (arr:set! (.velts v) 1 (* (arr:aref (.velts v) 1) a))
  (arr:set! (.velts v) 2 (* (arr:aref (.velts v) 2) a))
  v)

(declare vcross ((types:RuntimeRepr :t) (Num :t) => (Vec3 :t) -> (Vec3 :t) -> (Vec3 :t)))
(define (vcross a b)
  (let ((ax (vx a)) (ay (vy a)) (az (vz a))
        (bx (vx b)) (by (vy b)) (bz (vz b)))
    (v3 (- (* ay bz) (* by az))
        (- (* az bx) (* bz ax))
        (- (* ax by) (* bx ay)))))

(declare vadd! ((Num :t) => (Vec3 :t) -> (Vec3 :t) -> (Vec3 :t)))
(define (vadd! a b)
  (arr:set! (.velts a) 0 (+ (arr:aref (.velts a) 0) (arr:aref (.velts b) 0)))
  (arr:set! (.velts a) 1 (+ (arr:aref (.velts a) 1) (arr:aref (.velts b) 1)))
  (arr:set! (.velts a) 2 (+ (arr:aref (.velts a) 2) (arr:aref (.velts b) 2)))
  a)

(declare vsub! ((Num :t) => (Vec3 :t) -> (Vec3 :t) -> (Vec3 :t)))
(define (vsub! a b)
  (arr:set! (.velts a) 0 (- (arr:aref (.velts a) 0) (arr:aref (.velts b) 0)))
  (arr:set! (.velts a) 1 (- (arr:aref (.velts a) 1) (arr:aref (.velts b) 1)))
  (arr:set! (.velts a) 2 (- (arr:aref (.velts a) 2) (arr:aref (.velts b) 2)))
  a)

(declare vmul! ((Num :t) => (Vec3 :t) -> (Vec3 :t) -> (Vec3 :t)))
(define (vmul! a b)
  (arr:set! (.velts a) 0 (* (arr:aref (.velts a) 0) (arr:aref (.velts b) 0)))
  (arr:set! (.velts a) 1 (* (arr:aref (.velts a) 1) (arr:aref (.velts b) 1)))
  (arr:set! (.velts a) 2 (* (arr:aref (.velts a) 2) (arr:aref (.velts b) 2)))
  a)

(declare vdiv! ((Reciprocable :t) => (Vec3 :t) -> (Vec3 :t) -> (Vec3 :t)))
(define (vdiv! a b)
  (arr:set! (.velts a) 0 (/ (arr:aref (.velts a) 0) (arr:aref (.velts b) 0)))
  (arr:set! (.velts a) 1 (/ (arr:aref (.velts a) 1) (arr:aref (.velts b) 1)))
  (arr:set! (.velts a) 2 (/ (arr:aref (.velts a) 2) (arr:aref (.velts b) 2)))
  a)

(declare vnormalize! ((Radical :t) (Reciprocable :t) => (Vec3 :t) -> (Vec3 :t)))
(define (vnormalize! v)
  (vscale! v (reciprocal (sqrt (vdot v v)))))


(define-struct (Ray :t)
  (org "Origin"    (Vec3 :t))
  (dir "Direction" (Vec3 :t)))

(declare copy-ray (types:RuntimeRepr :t => (Ray :t) -> (Ray :t)))
(define (copy-ray r)
  (Ray (vdup (.org r)) (vdup (.dir r))))

;; map [0d0, 1d0] to [0, 255]
(declare clampu8 (Double-Float -> U8))
(define (clampu8 f)
  (let ((v (toInteger (floor (* f 255.5d0)))))
    (cond ((< v 0) 0)
          ((> v 255) 255)
          (true (fromInt v)))))

;; calculates basis vectors aligned to N
(declare orthoBasis ((Vec3 Double-Float) -> (Tuple3 (Vec3 Double-Float)
                                                    (Vec3 Double-Float)
                                                    (Vec3 Double-Float))))
(define (orthoBasis n)
  (let ((v (cond ((and (< -0.6d0 (vx n)) (< (vx n) 0.6d0)) (v3 1d0 0d0 0d0))
                 ((and (< -0.6d0 (vy n)) (< (vy n) 0.6d0)) (v3 0d0 1d0 0d0))
                 ((and (< -0.6d0 (vz n)) (< (vz n) 0.6d0)) (v3 0d0 0d0 1d0))
                 (true (v3 1d0 0d0 0d0))))
        (s (vnormalize! (vcross v n))))
    (Tuple3 s (vnormalize! (vcross n s)) n)))

;; returns a unit vector points to a random direction
(declare random-direction (Unit -> (Vec3 Double-Float)))
(define (random-direction)
  (let ((r (lisp Double-Float () (cl:random 1d0)))
        (phi (* 2.0d0 (* pi (lisp Double-Float () (cl:random 1d0)))))
        (rho (sqrt (- 1.0d0 r))))
    (v3 (* (cos phi) rho) (* (sin phi) rho) (sqrt r))))

;; Represents intersection.  One instance of intersection is allocated
;; per scanline and reused by all traces.
(define-struct Intersection
  (t           "distance from the ray origin"   (cell:Cell Double-Float))
  (updater     "a thunk to update intersection" (cell:Cell (Intersection
                                                            -> (Ray Double-Float)
                                                            -> Boolean)))
  (i-point     "where the ray intersects"       (Vec3 Double-Float))
  (i-normal    "surface normal"                 (Vec3 Double-Float))
  (col         "color of the surface"           (Vec3 Double-Float))
  (emissiveCol "emissive color of the surface"  (Vec3 Double-Float))
  (tmp         "temporary working area"         (Vec3 Double-Float))
  )

(declare update-noop (Intersection -> (Ray Double-Float) -> Boolean))
(define (update-noop _isect _r)
  false)

(declare make-intersection (Unit -> Intersection))
(define (make-intersection)
  (Intersection (cell:new 1d30)
                (cell:new update-noop)
                (v3 0d0 0d0 0d0)
                (v3 0d0 0d0 0d0)
                (v3 0d0 0d0 0d0)
                (v3 0d0 0d0 0d0)
                (v3 0d0 0d0 0d0)))

(declare reset-intersection! (Intersection -> Intersection))
(define (reset-intersection! isect)
  (cell:write! (.t isect) 1d30)
  (let _ = (cell:write! (.updater isect) update-noop))
  isect)

;;-------------------------------------------------
;; Primitives.
;; A primitive returns a procedure that calculates intersection.

(declare Sphere ((Vec3 Double-Float) -> Double-Float
                 -> (Vec3 Double-Float) -> (Vec3 Double-Float)
                 -> (Intersection -> (Ray Double-Float) -> Boolean)))
(define (Sphere center radius col emissiveCol)
  (let ((declare update-intersection! (Intersection -> (Ray Double-Float) -> Boolean))
        (update-intersection!
          (fn (isect ry)
            (let ((nn (.i-normal (the Intersection isect)))
                  (pp (.i-point (the Intersection isect))))
              (vcopy! pp (.dir ry))
              (vscale! pp (cell:read (.t isect)))
              (vadd! pp (.org ry))
              (vcopy! nn pp)
              (vsub! nn center)
              (vnormalize! nn)
              (vcopy! (.col isect) col)
              (vcopy! (.emissiveCol isect) emissiveCol)
              true)))

        (radius^2 (* radius radius))

        (declare intersect! (Intersection -> (Ray Double-Float) -> Boolean))
        (intersect!
          (fn (isect ry)
            (let ((v (.tmp (the Intersection isect))))
              (vcopy! v (.org ry))
              (vsub! v center)
              (let ((B (vdot v (.dir (the (Ray Double-Float) ry))))
                    (C (- (vdot v v) radius^2))
                    (D (- (* B B) C)))
                (if (> D 0d0)
                    (let ((t (- 0d0 (- B (sqrt D)))))
                      (if (and (< 0d0 t)
                               (< t (cell:read (.t isect))))
                          (progn
                            (cell:write! (.t isect) t)
                            (let _ = (cell:write! (.updater isect) update-intersection!))
                            true)
                          false))
                    false))))))
    intersect!))

(declare Plane ((Vec3 Double-Float) -> (Vec3 Double-Float)
                -> (Vec3 Double-Float) -> (Vec3 Double-Float)
                -> (Intersection -> (Ray Double-Float) -> Boolean)))
(define (Plane p n col emissiveCol)
  (let ((declare update-intersection! (Intersection -> (Ray Double-Float) -> Boolean))
        (update-intersection!
          (fn (isect ry)
            (vcopy! (.i-normal isect) n)
            (vcopy! (.i-point isect) (.dir (the (Ray Double-Float) ry)))
            (vscale! (.i-point isect) (cell:read (.t isect)))
            (vadd! (.i-point isect) (.org ry))
            (vcopy! (.col isect) col)
            (vcopy! (.emissiveCol isect) emissiveCol)
            true))

        (p.n (vdot p n))

        (declare intersect! (Intersection -> (Ray Double-Float) -> Boolean))
        (intersect!
          (fn (isect ry)
            (let ((v (vdot (.dir (the (Ray Double-Float) ry)) n)))
              (if (> (abs v) 1.0d-6)
                  (let ((t (/ (- p.n (vdot (.org ry) n)) v)))
                    (if (and (< 0d0 t)
                             (< t (cell:read (.t isect))))
                        (progn
                          (cell:write! (.t isect) t)
                          (let _ = (cell:write! (.updater isect) update-intersection!))
                          true)
                        false))
                  false))))
        )
    intersect!))

;; Find closest intersection of ray and objects.
;; Updates isect.  Returns true if ray intersects, false otherwise.
(declare find-is! (Intersection -> (Ray Double-Float)
                                -> (List (Intersection -> (Ray Double-Float) -> Boolean))
                                -> Unit))
(define (find-is! isect ray objs)
  (reset-intersection! isect)
  (let ((declare rec ((List (Intersection -> (Ray Double-Float) -> Boolean)) -> Unit))
        (rec
          (fn (objs)
            (match objs
              ((Cons obj rest)
               (obj isect ray)
               (rec rest))
              ((Nil) Unit)))))
    (rec objs)
    ((cell:read (.updater isect)) isect ray)
    Unit))

;;-------------------------------------------------
;; Main entry
;;

(define IMAGE_WIDTH (the UFix 256))
(define IMAGE_HEIGHT (the UFix 256))
(define NSUBSAMPLES (the UFix 2))
(define NPATH_SAMPLES (the UFix 128))
(define NAO_SAMPLES (the UFix 8))
(define MAX_TRACE_DEPTH (the UFix 16))

;; (define (run-ambient-occlusion)
;;   (do-scan ambient-occlusion))

(declare times (UFix -> (UFix -> Unit) -> Unit))
(define (times n proc)
  (let ((declare rec (UFix -> Unit))
        (rec (fn (i)
               (if (< i n)
                   (progn
                     (proc i)
                     (rec (+ i 1)))
                   Unit))))
    (rec 0)))

(declare do-scan (((Ray Double-Float) -> Intersection -> Double-Float) -> Unit))
(define (do-scan sampler)
  (times IMAGE_HEIGHT
         (fn (y)
           (render IMAGE_WIDTH IMAGE_HEIGHT y NSUBSAMPLES sampler)
           Unit)))

(declare render (UFix -> UFix -> UFix -> UFix
                      -> ((Ray Double-Float) -> Intersection -> Double-Float)
                      -> (arr:LispArray U8)))
(define (render width height y nsamp sampler)
  (let ((isect (make-intersection))
        (sum (cell:new 0d0))
        (img (lisp (arr:LispArray U8) (width)
               (cl:make-array (cl:list (cl:* width 4))
                              :element-type 'cl:fixnum)))
        (double (fn (n)
                  (lisp Double-Float (n) (cl:coerce n 'cl:double-float)))))

    (times
     width
     (fn (i)
       (cell:write! sum 0d0)
       (times
        nsamp
        (fn (u)
          (times
           nsamp
           (fn (v)
             (let ((px (/ (+ (double i)
                             (- (/ (double u) (double nsamp))
                                (/ (double width) 2d0)))
                          (/ (double width) 2d0)))
                   (py (- 0d0
                          (/ (+ (double y) (- (/ (double v) (double nsamp))
                                              (/ (double height) 2d0)))
                             (/ (double height) 2d0)))))
               (cell:write! sum
                            (+ (cell:read sum)
                               (sampler (Ray (v3 0d0 0d0 0d0)
                                             (vnormalize! (v3 px py -1d0)))
                                        isect)))
               Unit)))))
       (let ((area (lisp Double-Float (nsamp)
                     (cl:coerce (* nsamp nsamp) 'cl:double-float)))
             (val (cell:read sum)))
         (arr:set! img (+ (* i 4) 0) (clampu8 (/ val area)))
         (arr:set! img (+ (* i 4) 1) (clampu8 (/ val area)))
         (arr:set! img (+ (* i 4) 2) (clampu8 (/ val area))))))
    img))

#|

;;-------------------------------------------------
;; Ambient occlusion
;;

(define *objects-ao*
  (list
    (Sphere (vec -2.0 0.0 -3.5) 0.5 #f #f)
    (Sphere (vec -0.5 0.0 -3.0) 0.5 #f #f)
    (Sphere (vec 1.0 0.0 -2.2) 0.5  #f #f)
    (Plane (vec 0.0 -0.5 0.0) (vec 0.0 1.0 0.0) #f #f)))

(define (ambientOcclusion ray isect)
  (if (not (find-is! isect ray *objects-ao*))
    0.0
    (let ([eps 0.00001]
          [ntheta NAO_SAMPLES]
          [nphi NAO_SAMPLES])
      ;; set ray origin
      (let1 p (org ray)
        (vset! p (is-n isect)) (vmul! p eps) (vadd! p (is-p isect)))
      (let loop ([j 0] [i 0] [occ 0.0])
        (cond
         [(= j ntheta) (/ (- (* ntheta nphi) occ) (* ntheta nphi))]
         [(= i nphi)   (loop (+ j 1) 0 occ)]
         [else
          (receive (x y z) (random-direction)
            (receive (b0 b1 b2) (orthoBasis (is-n isect))
              ;; update ray direction.  be careful not to overwrite b2,
              ;; which is the same vector as (is-n isect).
              (let1 d (dir ray)
                (vmul! b0 x) (vmul! b1 y) (vset! d b2) (vmul! d z)
                (vadd! d b0) (vadd! d b1))
              (let1 hit (find-is! isect ray *objects-ao*)
                (loop j (+ i 1) (if hit (+ occ 1.0) occ)))))])))))


|#

)
