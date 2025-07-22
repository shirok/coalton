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

  (declare vmul-s! ((Num :t) => (Vec3 :t) -> :t -> (Vec3 :t)))
  (define (vmul-s! a s)
    (arr:set! (.velts a) 0 (* (arr:aref (.velts a) 0) s))
    (arr:set! (.velts a) 1 (* (arr:aref (.velts a) 1) s))
    (arr:set! (.velts a) 2 (* (arr:aref (.velts a) 2) s))
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
  (declare clampu8 (F64 -> U8))
  (define (clampu8 f)
    (let ((v (toInteger (floor (* f 255.5d0)))))
      (cond ((< v 0) 0)
            ((> v 255) 255)
            (true (fromInt v)))))

  ;; calculates basis vectors aligned to N
  (declare orthoBasis ((Vec3 F64) -> (Tuple3 (Vec3 F64)
                                             (Vec3 F64)
                                             (Vec3 F64))))
  (define (orthoBasis n)
    (let ((v (cond ((and (< -0.6d0 (vx n)) (< (vx n) 0.6d0)) (v3 1d0 0d0 0d0))
                   ((and (< -0.6d0 (vy n)) (< (vy n) 0.6d0)) (v3 0d0 1d0 0d0))
                   ((and (< -0.6d0 (vz n)) (< (vz n) 0.6d0)) (v3 0d0 0d0 1d0))
                   (true (v3 1d0 0d0 0d0))))
          (s (vnormalize! (vcross v n))))
      (Tuple3 s (vnormalize! (vcross n s)) n)))

  ;; returns a unit vector points to a random direction
  (declare random-direction (Unit -> (Vec3 F64)))
  (define (random-direction)
    (let ((r (lisp F64 () (cl:random 1d0)))
          (phi (* 2.0d0 (* pi (lisp F64 () (cl:random 1d0)))))
          (rho (sqrt (- 1.0d0 r))))
      (v3 (* (cos phi) rho) (* (sin phi) rho) (sqrt r))))

  ;; Represents intersection.  One instance of intersection is allocated
  ;; per scanline and reused by all traces.
  (define-struct Intersection
    (t           "distance from the ray origin"   (cell:Cell F64))
    (updater     "a thunk to update intersection" (cell:Cell (Intersection
                                                              -> (Ray F64)
                                                              -> Boolean)))
    (i-point     "where the ray intersects"       (Vec3 F64))
    (i-normal    "surface normal"                 (Vec3 F64))
    (col         "color of the surface"           (Vec3 F64))
    (emissiveCol "emissive color of the surface"  (Vec3 F64))
    (tmp         "temporary working area"         (Vec3 F64))
    )

  (declare update-noop (Intersection -> (Ray F64) -> Boolean))
  (define (update-noop _isect _r)
    False)

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

  (declare Sphere ((Vec3 F64) -> F64
                   -> (Vec3 F64) -> (Vec3 F64)
                   -> (Intersection -> (Ray F64) -> Boolean)))
  (define (Sphere center radius col emissiveCol)
    (let ((declare update-intersection! (Intersection -> (Ray F64) -> Boolean))
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

          (declare intersect! (Intersection -> (Ray F64) -> Boolean))
          (intersect!
            (fn (isect ry)
              (let ((v (.tmp (the Intersection isect))))
                (vcopy! v (.org ry))
                (vsub! v center)
                (let ((B (vdot v (.dir (the (Ray F64) ry))))
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

  (declare Plane ((Vec3 F64) -> (Vec3 F64)
                  -> (Vec3 F64) -> (Vec3 F64)
                  -> (Intersection -> (Ray F64) -> Boolean)))
  (define (Plane p n col emissiveCol)
    (let ((declare update-intersection! (Intersection -> (Ray F64) -> Boolean))
          (update-intersection!
            (fn (isect ry)
              (vcopy! (.i-normal isect) n)
              (vcopy! (.i-point isect) (.dir (the (Ray F64) ry)))
              (vscale! (.i-point isect) (cell:read (.t isect)))
              (vadd! (.i-point isect) (.org ry))
              (vcopy! (.col isect) col)
              (vcopy! (.emissiveCol isect) emissiveCol)
              true))

          (p.n (vdot p n))

          (declare intersect! (Intersection -> (Ray F64) -> Boolean))
          (intersect!
            (fn (isect ry)
              (let ((v (vdot (.dir (the (Ray F64) ry)) n)))
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
  (declare find-is! (Intersection -> (Ray F64)
                                  -> (List (Intersection -> (Ray F64) -> Boolean))
                                  -> Boolean))
  (define (find-is! isect ray objs)
    (reset-intersection! isect)
    (let ((declare %rec ((List (Intersection -> (Ray F64) -> Boolean)) -> Unit))
          (%rec
            (fn (objs)
              (match objs
                ((Cons obj rest)
                 (obj isect ray)
                 (%rec rest))
                ((Nil) Unit)))))
      (%rec objs)
      ((cell:read (.updater isect)) isect ray)))

  ;;-------------------------------------------------
  ;; Main entry
  ;;

  (define (run-ambient-occlusion)
    (do-scan ambientOcclusion))

  (define IMAGE_WIDTH (the UFix 256))
  (define IMAGE_HEIGHT (the UFix 256))
  (define NSUBSAMPLES (the UFix 2))
  (define NPATH_SAMPLES (the UFix 128))
  (define NAO_SAMPLES (the UFix 8))
  (define MAX_TRACE_DEPTH (the UFix 16))

  (declare times (UFix -> (UFix -> Unit) -> Unit))
  (define (times n proc)
    (let ((declare %rec (UFix -> Unit))
          (%rec (fn (i)
                  (if (< i n)
                      (progn
                        (proc i)
                        (%rec (+ i 1)))
                      Unit))))
      (%rec 0)))

  (declare do-scan (((Ray F64) -> Intersection -> F64) -> Unit))
  (define (do-scan sampler)
    (l:dotimes (y IMAGE_HEIGHT)
      (let ((scanline (render IMAGE_WIDTH IMAGE_HEIGHT y NSUBSAMPLES sampler)))
        (lisp Unit (y)
          (cl:format cl:t "LINE ~d~%" y)
          Unit))))

  (declare render (UFix -> UFix -> UFix -> UFix
                        -> ((Ray F64) -> Intersection -> F64)
                        -> (arr:LispArray U8)))
  (define (render width height y nsamp sampler)
    (let ((isect (make-intersection))
          (sum (cell:new 0d0))
          (img (lisp (arr:LispArray U8) (width)
                 (cl:make-array (cl:list (cl:* width 4))
                                :element-type '(cl:unsigned-byte 8))))
          (double (fn (n) (the F64 (fromInt (into n))))))

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
         (let ((area (double (* nsamp nsamp)))
               (val (cell:read sum)))
           (arr:set! img (+ (* i 4) 0) (clampu8 (/ val area)))
           (arr:set! img (+ (* i 4) 1) (clampu8 (/ val area)))
           (arr:set! img (+ (* i 4) 2) (clampu8 (/ val area))))))
      img))

  ;; ambient occlusion scanner
  (define *objects-ao*
    (make-list
     (Sphere (v3 -2.0d0 0.0d0 -3.5d0) 0.5d0
             (v3 1.0d0 1.0d0 1.0d0) (v3 0.0d0 0.0d0 0.0d0))
     (Sphere (v3 -0.5d0 0.0d0 -3.0d0) 0.5d0
             (v3 1.0d0 1.0d0 1.0d0) (v3 0.0d0 0.0d0 0.0d0))
     (Sphere (v3 1.0d0 0.0d0 -2.2d0) 0.5d0
             (v3 1.0d0 1.0d0 1.0d0) (v3 0.0d0 0.0d0 0.0d0))
     (Plane (v3 0.0d0 -0.5d0 0.0d0) (v3 0.0d0 1.0d0 0.0d0)
            (v3 1.0d0 1.0d0 1.0d0) (v3 0.0d0 0.0d0 0.0d0))))

  (define (ambientOcclusion ray isect)
    (if (not (find-is! isect ray *objects-ao*))
        0.0d0
        (let ((eps 0.00001d0)
              (ntheta NAO_SAMPLES)
              (nphi NAO_SAMPLES)
              (p (.org ray)))
          ;; set ray origin
          (vcopy! p (.i-normal isect))
          (vmul-s! p eps)
          (vadd! p (.i-point isect))
          (rec %loop ((j 0) (i 0) (occ 0.0d0))
            (cond
              ((== j ntheta) (/ (- (the F64 (fromInt (into (* ntheta nphi)))) occ)
                                (the F64 (fromInt (into (* ntheta nphi))))))
              ((== i nphi)   (%loop (+ j 1) 0 occ))
              (True
               (let dir = (random-direction))
               (let (Tuple3 b0 b1 b2) = (orthoBasis (.i-normal isect)))
               ;; update ray direction.  be careful not to overwrite b2,
               ;; which is the same vector as (.i-normal isect).
               (let d = (.dir ray))
               (vmul-s! b0 (vx dir))
               (vmul-s! b1 (vy dir))
               (vcopy! d b2)
               (vmul-s! d (vz dir))
               (vadd! d b0)
               (vadd! d b1)
               (let hit = (find-is! isect ray *objects-ao*))
               (%loop j (+ i 1) (if hit (+ occ 1.0d0) occ))))))))


  ;; path tracing scanner
#|
  (define *objects-pt*
    (make-list
     (Sphere (v3 -1.05d0 0.0d0 -2.0d0) 0.5d0
             (v3 0.75d0 0.0d0 0.0d0) (v3 0.0d0 0.0d0 0.0d0))
     (Sphere (v3 0.0d0 0.0d0 -2.0d0) 0.5d0
             (v3 1.0d0 1.0d0 1.0d0) (v3 1.0d0 1.0d0 1.0d0))
     (Sphere (v3 1.05d0 0.0d0 -2.0d0) 0.5d0
             (v3 0.0d0 0.0d0 1.0d0) (v3 0.0d0 0.0d0 0.0d0))
     (Plane (v3 0.0d0 -0.5d0 0.0d0) (v3 0.0d0 1.0d0 0.0d0)
            (v3 1.0d0 1.0d0 1.0d0) (v3 0.0d0 0.0d0 0.0d0))))

  (declare trace! ((Ray F64) -> Intersection -> F64 -> (Vec3 F64)))
  (define (trace! r isect depth)
    (cond ((not (find-is! isect r *objects-pt*))
           (v3 0.7d0 0.7d0 0.7d0))
          (True
           (let dir = (random-direction))
           (let (Tuple3 b0 b1 b2) = (orthoBasis (.i-normal isect)))
           ;; update ray direction.  be careful not to overwrite b2,
           ;; which is the same vector as (.i-normal isect).
           (vmul-s! b0 (vx dir))
           (vmul-s! b1 (vy dir))
           (vcopy! (.dir r) b2)
           (vmul-s! (.dir r) (vz dir))
           (vadd! (.dir r) b0)
           (vadd! (.dir r) b1)
           ;; update ray origin
           (vcopy! (.org r) (.dir r))
           (vmul-s! (.org r) 0.00001d0)
           (vadd! (.org r) (.i-point isect))
           (let col = (if (== depth (- MAX_TRACE_DEPTH 1))
                          (v3 0.0d0 0.0d0 0.0d0)
                          (trace! r isect (+ depth 1))))
           (vmul! col (.col isect))
           (vadd! col (.emissiveCol isect))
           col)))

  (declare pathTrace ((Ray F64) -> Intersection -> (Vec3 F64)))
  (define (pathTrace r isect)
    (let col = (v3 0.0d0 0.0d0 0.0d0))
    (dotimes (i NPATH_SAMPLES)
      (vadd! col (trace! (copy-ray r) isect 0)))
    (vmul! col (/ 1.0d0 NPATH_SAMPLES))
    col)
|#
  )
