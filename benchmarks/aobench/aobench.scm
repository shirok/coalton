;; based on http://lucille.atso-net.jp/aobench/
;; and http://d.hatena.ne.jp/mjt/20090206/p1
;;     http://d.hatena.ne.jp/mjt/20090209/p1

(use gauche.uvector)
(use gauche.threads)
(use gauche.parseopt)
(use srfi-1)
(use srfi-27)
(use srfi-42)
(use math.const)
(use util.match)

;;-------------------------------------------------
;; Common utilities
;;

;; vector stuff
(define vec f64vector)
(define-inline (vx v) (f64vector-ref v 0))
(define-inline (vy v) (f64vector-ref v 1))
(define-inline (vz v) (f64vector-ref v 2))

(define-macro (vset! v . args)
  (match args
    [(x y z) `(begin (f64vector-set! ,v 0 ,x)
                     (f64vector-set! ,v 1 ,y)
                     (f64vector-set! ,v 2 ,z))]
    [(vv)    `(f64vector-copy! ,v ,vv)]))

(define vdot f64vector-dot)
(define vscale f64vector-mul)
(define (vcross a b)
  (let ([ax (vx a)] [ay (vy a)] [az (vz a)]
        [bx (vx b)] [by (vy b)] [bz (vz b)])
    (vec (- (* ay bz) (* by az))
         (- (* az bx) (* bz ax))
         (- (* ax by) (* bx ay)))))

(define vadd! f64vector-add!)
(define vsub! f64vector-sub!)
(define vmul! f64vector-mul!)
(define vdiv! f64vector-div!)
(define (vnormalize! a)
  (let1 len (sqrt (vdot a a))
    (when (> len 1.0e-6) (vdiv! a len)))
  a)

;; ray
(define org car)
(define dir cdr)
(define new-ray cons)
(define (copy-ray ray)
  (new-ray (f64vector-copy (org ray)) (f64vector-copy (dir ray))))

;; clamp for final pixel value
(define (clampu8 f)
  (clamp (floor->exact (* f 255.5)) 0 255))

;; calculates basis vectors aligned to N
(define (orthoBasis n)
  (let* ([v (cond [(< -0.6 (vx n) 0.6) '#f64(1.0 0.0 0.0)]
                  [(< -0.6 (vy n) 0.6) '#f64(0.0 1.0 0.0)]
                  [(< -0.6 (vz n) 0.6) '#f64(0.0 0.0 1.0)]
                  [else '#f64(vec 1.0 0.0 0.0)])]
         [s (vnormalize! (vcross v n))])
    (values s (vnormalize! (vcross n s)) n)))

;; returns a unit vector points to a random direction
(define (random-direction)
  (let ([r (random-real)] [phi (* 2.0 pi (random-real))])
    (let1 rho (sqrt (- 1.0 r))
      (values (* (cos phi) rho) (* (sin phi) rho) (sqrt r)))))

;; Macro define-simple-struct creates a bunch of functions and macros
;; to emulate a structure by a vector.
;; (This is a kludge.  Future version of Gauche will have records as
;; efficient as this trick.)
;;
;; NAME is a symbol to name the structure type.  TAG is some value
;; (usually a symbol or an integer) to indicate the type of the
;; structure.
;;
;; (define-simple-struct <name> <tag> <constructor> [(<slot-spec>*)])
;;
;; <constructor> : <symbol> | #f
;; <slot-spec>   : <slot-name> | (<slot-name> [<init-value>])
;;
;; For each <slot-spec>, the following accessor/modifier are automatially
;; generated.  
;;
;;   NAME-SLOT      - accessor (macro)
;;   NAME-SLOT-set! - modifier (macro)

(define-macro (define-simple-struct name tag constructor . opts)
  (let-optionals* opts ((slot-defs '()))
    (define (make-constructor)
      (let ((args (gensym))
            (num-slots  (length slot-defs))
            (slot-names (map (lambda (s) (if (symbol? s) s (car s)))
                             slot-defs))
            (init-vals  (map (lambda (s) (if (symbol? s) #f (cadr s)))
                             slot-defs)))
        `(define-macro (,constructor . ,args)
           (match ,args
             ,@(let loop ((n 0)
                          (r '()))
                 (if (> n num-slots)
                   r
                   (let1 carg (take slot-names n)
                     (loop (+ n 1)
                           (cons
                            `(,carg
                              (list 'vector
                                    ,@(if tag `(',tag) '())
                                    ,@carg
                                    ,@(map (cut list 'quote <>)
                                           (list-tail init-vals n))))
                            r)))
                   ))))
        ))
    `(begin
       ,@(if constructor
           `(,(make-constructor))
           '())
       ,@(let loop ((s slot-defs) (i (if tag 1 0)) (r '()))
           (if (null? s)
             (reverse r)
             (let* ((slot-name (if (pair? (car s)) (caar s) (car s)))
                    (acc (string->symbol #`",|name|-,|slot-name|"))
                    (mod (string->symbol #`",|name|-,|slot-name|-set!")))
               (loop (cdr s)
                     (+ i 1)
                     (list*
                      `(define-macro (,acc obj)
                         `(vector-ref ,obj ,,i))
                      `(define-macro (,mod obj val)
                         `(vector-set! ,obj ,,i ,val))
                      r))))))
    ))

;; Represents intersection.  One instance of intersection is allocated
;; per scanline and reused by all traces.
(define-simple-struct is 'intersection make-intersection
  ([t   1.0e+30]                ; distance from the ray origin
   [updater #f]                 ; a thunk to update intersection
   [p   (make-f64vector 3 0.0)] ; where the ray intersects
   [n   (make-f64vector 3 0.0)] ; surface normal
   [col #f]                     ; color of the surface (pathtrace)
   [emissiveCol #f]             ; emissive color of the surface (pathtrace)
   [tmp (make-f64vector 3 0.0)] ; temporary working area
                                ;  (ugly, but to reduce allocation)
   ))

(define (reset-intersection! is)
  (is-t-set! is 1.0e+30)
  (is-updater-set! is #f)
  is)

;;-------------------------------------------------
;; Primitives.
;; A primitive returns a procedure that calculates intersection.

(define (Sphere center radius col emissiveCol)

  (define (update-intersection! isect ray)
    (let ([pp (is-p isect)]
          [nn (is-n isect)])
      (vset! pp (dir ray))
      (vmul! pp (is-t isect))
      (vadd! pp (org ray))
      (vset! nn pp)
      (vsub! nn center)
      (vnormalize! nn))
    (is-col-set! isect col)
    (is-emissiveCol-set! isect emissiveCol))

  (define radius^2 (* radius radius))
  
  (lambda (isect ray)
    (let* ([rs (rlet1 v (is-tmp isect)
                 (vset! v (org ray)) (vsub! v center))]
           [B (vdot rs (dir ray))]
           [C (- (vdot rs rs) radius^2)]
           [D (- (* B B) C)])
      (when (> D 0.0)
        (let1 t (- 0.0 B (sqrt D))
          (when (< 0.0 t (is-t isect))
            (is-t-set! isect t)
            (is-updater-set! isect update-intersection!)))))))

(define (Plane p n col emissiveCol)

  (define (update-intersection! isect ray)
    (vset! (is-n isect) n)
    (vset! (is-p isect) (dir ray))
    (vmul! (is-p isect) (is-t isect))
    (vadd! (is-p isect) (org ray))
    (is-col-set! isect col)
    (is-emissiveCol-set! isect emissiveCol))

  (define p.n (vdot p n))

  (lambda (isect ray)
    (let1 v (vdot (dir ray) n)
      (when (> (abs v) 1.0e-6)
	(let1 t (/ (- (- (vdot (org ray) n) p.n)) v)
	  (when (< 0.0 t (is-t isect))
	    (is-t-set! isect t)
	    (is-updater-set! isect update-intersection!)))))))

;; Find closest intersection of ray and objects.
;; Updates isect.  Returns true if ray intersects, false otherwise.
(define (find-is! isect ray objs)
  (reset-intersection! isect)
  (let loop ((objs objs))
    (unless (null? objs)
      ((car objs) isect ray)
      (loop (cdr objs))))
  (cond [(is-updater isect) => (lambda (u) (u isect ray) #t)]
        [else #f]))

;;-------------------------------------------------
;; Main entry
;;

(define-constant IMAGE_WIDTH 256)
(define-constant IMAGE_HEIGHT 256)
(define-constant NSUBSAMPLES 2)
(define-constant NPATH_SAMPLES 128)
(define-constant NAO_SAMPLES 8)
(define-constant MAX_TRACE_DEPTH 16)

(define (main args)
  (let-args (cdr args) ([type   "t=y" 'ao]
                        [nprocs "p=i" 1]
                        [else => (lambda _ (usage))])
    (call-with-output-file "out.raw"
      (cute (if (> nprocs 1) do-parallel do-serial)
            (case type
              [(ao) ambientOcclusion]
              [(pt) pathTrace]
              [else (exit 1 "invalid type (must be either ao or pt): ~a" type)])
            nprocs <>))
    (print "DONE")
    0))

(define (usage)
  (print "Usage: aobench [-t ao|pt] [-p nprocs]")
  (exit 0))

(define (do-serial sampler _ out)
  (dotimes [y IMAGE_HEIGHT]
    (let1 scanline (render IMAGE_WIDTH IMAGE_HEIGHT y NSUBSAMPLES sampler)
      (print "LINE"y" - " (u8vector-length scanline))
      (write-block scanline out))))

;; naive parallelization by static partitioning
(define (do-parallel sampler nproc out)
  (let1 lines (make-vector IMAGE_HEIGHT #f)
    (define (calc y)
      (print "LINE"y)
      (set! (ref lines y)
            (render IMAGE_WIDTH IMAGE_HEIGHT y NSUBSAMPLES sampler)))
    (define (task off)
      (do-ec (: y off IMAGE_HEIGHT nproc) (calc y)))
    (for-each thread-join!
              (list-ec (: off nproc)
                       (thread-start! (make-thread (cut task off)))))
    (dotimes [y IMAGE_HEIGHT] (write-block (ref lines y) out))))

(define (render width height y nsamp sampler)
  (define isect (make-intersection))
  (define sum (vec 0.0 0.0 0.0))
  
  (rlet1 img (make-u8vector (* width 4))
    (dotimes [i width]
      (vset! sum 0 0 0)
      (dotimes [u nsamp]
        (dotimes [v nsamp]
          (let ([px (/ (+ i (- (/. u nsamp) (/. width 2))) (/. width 2))]
                [py (- (/ (+ y (- (/. v nsamp) (/. height 2))) (/. height 2)))])
            (vadd! sum (sampler (new-ray (vec 0.0 0.0 0.0)
                                         (vnormalize! (vec px py -1.0)))
                                isect)))))
      (u8vector-set! img (+ (* i 4) 0) (clampu8 (/ (vx sum) (* nsamp nsamp))))
      (u8vector-set! img (+ (* i 4) 1) (clampu8 (/ (vy sum) (* nsamp nsamp))))
      (u8vector-set! img (+ (* i 4) 2) (clampu8 (/ (vz sum) (* nsamp nsamp)))))))

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

;;-------------------------------------------------
;; Path tracer
;;

(define *objects-pt*
  (list
   (Sphere (vec -1.05 0.0 -2.0) 0.5 (vec 0.75 0.0 0.0) (vec 0.0 0.0 0.0))
   (Sphere (vec 0.0 0.0 -2.0) 0.5 (vec 1.0 1.0 1.0) (vec 1.0 1.0 1.0))
   (Sphere (vec 1.05 0.0 -2.0) 0.5 (vec 0.0 0.0 1.0) (vec 0.0 0.0 0.0))
   (Plane (vec 0.0 -0.5 0.0) (vec 0.0 1.0 0.0) (vec 1.0 1.0 1.0) (vec 0.0 0.0 0.0))))

(define (trace! ray isect depth)
  (if (not (find-is! isect ray *objects-pt*))
    (vec 0.7 0.7 0.7)
    (receive (x y z) (random-direction)
      (receive (b0 b1 b2) (orthoBasis (is-n isect))
        (let ([d (dir ray)] [o (org ray)])
          ;; update ray direction.  be careful not to overwrite b2,
          ;; which is the same vector as (is-n isect).
          (vmul! b0 x) (vmul! b1 y) (vset! d b2) (vmul! d z)
          (vadd! d b0) (vadd! d b1)
          ;; update ray origin
          (vset! o d) (vmul! o 0.00001) (vadd! o (is-p isect)))
        (let ([fr (is-col isect)]
              [Le (is-emissiveCol isect)])
          (rlet1 col (if (= depth (- MAX_TRACE_DEPTH 1))
                       (vec 0.0 0.0 0.0)
                       (trace! ray isect (+ depth 1)))
            (vmul! col fr) (vadd! col Le)))))))

(define (pathTrace ray isect)
  (rlet1 col (vec 0.0 0.0 0.0)
    (dotimes [i NPATH_SAMPLES]
      (vadd! col (trace! (copy-ray ray) isect 0)))
    (vmul! col (/. NPATH_SAMPLES))))

