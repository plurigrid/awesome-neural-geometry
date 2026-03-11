;;; world:// hamming a-z — path invariance + fibonacci automorphism
;;; parentheses coloring (Hy / Hylang)
;;;
;;; Finds the first sufficient automorphism of the Hamming distance
;;; matrix that bootstraps a consistent parentheses coloring scheme,
;;; using fibonacci-indexed letter subsets as seeds.

(import itertools)

;; --- encoding ---

(defn to-digits [val base width]
  (setv result [])
  (setv v val)
  (for [_ (range width)]
    (.insert result 0 (% v base))
    (setv v (// v base)))
  (tuple result))

(defn hamming [a b]
  (sum (gfor [x y] (zip a b) (if (!= x y) 1 0))))

;; --- letter setup ---

(setv LETTERS (lfor i (range 26) (chr (+ (ord "a") i))))

(defn encode-all [base width]
  (dfor ch LETTERS ch (to-digits (- (ord ch) (ord "a")) base width)))

(setv BIN-CODES (encode-all 2 5))
(setv TRI-CODES (encode-all 3 3))

(defn dist-matrix [codes]
  (setv n (len LETTERS))
  (setv mat (lfor _ (range n) (lfor _ (range n) 0)))
  (for [i (range n)]
    (for [j (range n)]
      (setv (get (get mat i) j)
            (hamming (get codes (get LETTERS i))
                     (get codes (get LETTERS j))))))
  mat)

;; --- fibonacci subset ---

(defn fib-indices [limit]
  (setv fibs [])
  (setv a 1)
  (setv b 1)
  (while (< a limit)
    (when (not-in a fibs)
      (.append fibs a))
    (setv tmp b)
    (setv b (+ a b))
    (setv a tmp))
  fibs)

;; --- automorphism detection ---
;; Swap automorphisms (transpositions preserving distance matrix)

(defn is-automorphism? [perm mat]
  (setv n (len perm))
  (for [i (range n)]
    (for [j (range (+ i 1) n)]
      (when (!= (get (get mat i) j)
                (get (get mat (get perm i)) (get perm j)))
        (return False))))
  True)

(defn find-swap-automorphisms [mat n]
  (setv autos [])
  (for [i (range n)]
    (for [j (range (+ i 1) n)]
      (setv perm (list (range n)))
      (setv (get perm i) j)
      (setv (get perm j) i)
      (when (is-automorphism? perm mat)
        (.append autos [i j]))))
  autos)

;; --- row-sum equivalence: letters with identical distance profiles ---
;; Two letters are "Hamming-equivalent" if they have the same row sum
;; (same total distance to all others). This is the first sufficient
;; automorphism class for bootstrapping paren coloring.

(defn row-sum-orbits [mat n]
  (setv by-sum {})
  (for [i (range n)]
    (setv s (sum (get mat i)))
    (when (not-in s by-sum)
      (setv (get by-sum s) []))
    (.append (get by-sum s) i))
  by-sum)

;; --- full distance-profile equivalence ---
;; Even stronger: two letters are profile-equivalent if their entire
;; sorted distance row is identical (not just the sum).

(defn profile-orbits [mat n]
  (setv by-profile {})
  (for [i (range n)]
    (setv prof (tuple (sorted (get mat i))))
    (when (not-in prof by-profile)
      (setv (get by-profile prof) []))
    (.append (get by-profile prof) i))
  by-profile)

;; --- parentheses coloring via orbits ---

(defn build-orbits [swaps n]
  (setv parent (list (range n)))
  (defn find-root [x]
    (while (!= (get parent x) x)
      (setv (get parent x) (get parent (get parent x)))
      (setv x (get parent x)))
    x)
  (for [[i j] swaps]
    (setv ri (find-root i))
    (setv rj (find-root j))
    (when (!= ri rj)
      (setv (get parent ri) rj)))
  (setv orbits {})
  (for [i (range n)]
    (setv r (find-root i))
    (when (not-in r orbits)
      (setv (get orbits r) []))
    (.append (get orbits r) i))
  orbits)

(setv PAREN-STYLES ["()" "[]" "{}" "<>" "«»"])

(defn color-parens [orbits]
  (setv coloring {})
  (for [[idx [root members]] (enumerate (.items orbits))]
    (setv style (get PAREN-STYLES (% idx (len PAREN-STYLES))))
    (for [m members]
      (setv (get coloring m) style)))
  coloring)

;; --- path cost ---

(defn path-cost [mat path-indices]
  (if (< (len path-indices) 2)
    0
    (sum (gfor i (range (- (len path-indices) 1))
           (get (get mat (get path-indices i)) (get path-indices (+ i 1)))))))

;; --- separator ---

(defn sep []
  (print (* "=" 70)))

;; --- main ---

(defn main []
  (sep)
  (print "  world:// hamming a-z — PATH INVARIANCE + PAREN COLORING (Hy)")
  (sep)

  (for [[label codes base width] [["Binary(5-bit)" BIN-CODES 2 5]
                                   ["TrueTrits(3-trit)" TRI-CODES 3 3]]]
    (setv mat (dist-matrix codes))
    (setv n 26)

    (print)
    (sep)
    (print (.format "  {}: Automorphisms & Parentheses Coloring" label))
    (sep)

    ;; find swap automorphisms
    (setv swaps (find-swap-automorphisms mat n))
    (print (.format "\n  Swap automorphisms (transpositions): {}" (len swaps)))

    ;; row-sum equivalence classes
    (setv rs-orbits (row-sum-orbits mat n))
    (print (.format "\n  Row-sum equivalence classes: {}" (len rs-orbits)))
    (for [[s members] (sorted (.items rs-orbits))]
      (setv letters-in (.join "" (lfor m members (get LETTERS m))))
      (print (.format "    sum={:3d}: {}  (size {})" s letters-in (len members))))

    ;; profile equivalence classes (first sufficient automorphism)
    (setv pr-orbits (profile-orbits mat n))
    (setv pr-list (sorted (.values pr-orbits) :key (fn [v] (get v 0))))
    (print (.format "\n  Distance-profile equivalence classes: {} (first sufficient automorphism)"
                    (len pr-orbits)))

    ;; use profile orbits for coloring (the finest meaningful partition)
    (setv orbits {})
    (for [[idx members] (enumerate pr-list)]
      (setv (get orbits idx) members))
    (setv coloring (color-parens orbits))
    (setv n-orbits (len orbits))
    (for [[idx members] (sorted (.items orbits))]
      (setv letters-in (.join "" (lfor m members (get LETTERS m))))
      (setv style (get coloring (get members 0)))
      (print (.format "    orbit {:2d}: {:12s}  paren: {}  size={}"
                      idx letters-in style (len members))))

    ;; colored alphabet
    (print "\n  Colored alphabet (profile-equivalence coloring):")
    (setv colored-parts
      (lfor i (range n)
        (+ (get (get coloring i) 0) (get LETTERS i) (get (get coloring i) 1))))
    (print (.format "    {}" (.join " " colored-parts)))

    ;; fibonacci subset path invariance
    (setv fib-idx (fib-indices 26))
    (setv fib-letters (lfor i fib-idx (get LETTERS i)))
    (print (.format "\n  Fibonacci subset: {}  indices={}"
                    (.join "" fib-letters) fib-idx))

    (setv fwd-cost (path-cost mat fib-idx))
    (setv rev-cost (path-cost mat (list (reversed fib-idx))))
    (print (.format "    Forward path cost:  {}" fwd-cost))
    (print (.format "    Reverse path cost:  {}" rev-cost))
    (print (.format "    Palindrome-invariant: {}" (= fwd-cost rev-cost)))

    ;; fibonacci-seeded coloring
    (print "\n  Fibonacci-seeded paren coloring:")
    (setv fib-set (set fib-idx))
    (setv fib-parts
      (lfor i (range n)
        (+ (if (in i fib-set) "*" " ")
           (get (get coloring i) 0)
           (get LETTERS i)
           (get (get coloring i) 1))))
    (print (.format "    {}" (.join " " fib-parts)))
    (print "    (* = fibonacci-indexed letters)")

    ;; full path comparison
    (print "\n  Path costs (natural a->z vs reverse z->a):")
    (setv nat-idx (list (range 26)))
    (setv rev-idx (list (reversed nat-idx)))
    (setv nat-cost (path-cost mat nat-idx))
    (setv rev-cost-full (path-cost mat rev-idx))
    (print (.format "    Natural (a->z):  {}" nat-cost))
    (print (.format "    Reverse (z->a):  {}" rev-cost-full))
    (print (.format "    Invariant: {}" (= nat-cost rev-cost-full)))

    ;; per-orbit path costs
    (print "\n  Per-orbit path costs:")
    (for [[root members] (sorted (.items orbits))]
      (when (> (len members) 1)
        (setv style (get coloring (get members 0)))
        (setv letters-in (.join "" (lfor m members (get LETTERS m))))
        (setv fwd (path-cost mat members))
        (setv rev (path-cost mat (list (reversed members))))
        (print (.format "    {} orbit [{}]:  fwd={}  rev={}  invariant={}"
                        style letters-in fwd rev (= fwd rev))))))

  ;; comparison
  (print)
  (sep)
  (print "  COMPARISON SUMMARY")
  (sep)

  (setv bin-mat (dist-matrix BIN-CODES))
  (setv tri-mat (dist-matrix TRI-CODES))
  (setv bin-pr (profile-orbits bin-mat 26))
  (setv tri-pr (profile-orbits tri-mat 26))
  (setv bin-rs (row-sum-orbits bin-mat 26))
  (setv tri-rs (row-sum-orbits tri-mat 26))

  (print (.format "\n  Binary profile-equiv classes:  {}" (len bin-pr)))
  (print (.format "  Trits  profile-equiv classes:  {}" (len tri-pr)))
  (print (.format "  Binary row-sum classes:         {}" (len bin-rs)))
  (print (.format "  Trits  row-sum classes:         {}" (len tri-rs)))
  (print "\n  Fewer classes = more symmetry = fewer paren colors needed")
  (print (.format "  Binary needs {} colors, Trits needs {} colors"
                  (len bin-pr) (len tri-pr)))
  (print "\n  Binary path a->z: 47  Trits path a->z: 35")
  (print "  Both are palindrome-invariant (forward == reverse)")
  (print "  Fibonacci path: binary=12, trits=10 (both palindrome-invariant)"))

(main)
