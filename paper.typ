#import "template.typ": *
#show: template.with(font-mode: "computer")

#import "@preview/hydra:0.6.3": hydra
#let paper-title = [Convergence analysis of teacher-model OBS methods in LR]
#set page(header: context {
  set text(size: 0.9em)
  emph(paper-title)
  h(1fr)
  emph(hydra(2))
})

#import "@preview/theoretic:0.4.0"
#import theoretic.presets.basic: *
#show ref: theoretic.show-ref

#let proof = proof.with(suffix: sym.square)

#let lemma(..args) = tang(theoretic.presets.basic.lemma(..args))
#let proposition(..args) = tang(theoretic.presets.basic.proposition(..args))
#let theorem(..args) = tang(theoretic.presets.basic.theorem(..args))

#let private(it) = text(gray, it)

#set text(size: 10pt)
#set heading(numbering: "1.1.1.")

#show heading.where(level:1): it => {
  counter(math.equation).update(0)
  it
}

// Do section-based equation numbering
#let num-eq(it) = {
  set math.equation(numbering: it => {
    let count = counter(heading.where(level: 1)).at(here()).first()
    numbering("(1.1)", count, it)
  })
  it
}

// Don't add "Equation" supplement before equation references
#show ref: it => {
  let eq = math.equation
  let el = it.element
  // Skip all other references.
  if el == none or el.func() != eq { return it }
  // Override equation references.
  link(el.location(), counter(eq).display(at: el.location()))
}

#align(center)[
  #text(size: 1.1em)[#strong(paper-title)]
  
  Matthew Ward and Kevin Miller
]

#to-do[Abstract]
= Introduction <sec-introduction>

#to-do[
  - Kaczmarz method
  - Brief introduction to online sample (batch?) selection
  - Justify reasoning behind teacher-informed methods
    - Real situations in which we might have a teacher available
    - Derivation of optimal criterion given perfect teacher
    - Define teacher-optimal criterion
  - Summarize results of paper, especially computational results  
]

== Setup
#let argmax = math.op("argmax", limits: true)

#let bx = $bold(x)$
#let ba = $bold(a)$
#let bb = $bold(b)$
#let beps = $bold(epsilon)$
With some $n, d in NN$ where $n > d$, we define the overdetermined linear system 
#num-eq[$
    A bx^* = bb + beps^*, wide A = vec(ba_1^trp, dots.v, ba_n^trp) in RR^(n times d), wide bx^* in RR^d, wide bb, beps^* in RR^n,
$ <eq:orig-problem>]

and we assume that the rows of $A$ are normalized for notational convenience: 
#num-eq[$
     forall i in {1, ..., n}, quad norm(ba_i^trp)_2 = 1.
$ <eq:normalized-rows>]

(Note that throughout this work, we use the asterisk ($*$) to denote optimality, not the matrix transpose.)

#let cardC = $abs(cal(C))$
We also assume that the noise vector $beps^*$ is sufficiently sparse to ensure that the noiseless rows of the system alone would guarantee a unique solution. That is, fix a "clean" set $cal(C) subset {1, .., n}$ with $cardC > d$, and assume $beps^*_i = 0$ for all $i in cal(C)$. Because $cardC > d$, we may assume that the matrix $A$ restricted to its clean rows, which we call $A_cal(C)$, has full column rank.

== Kaczmarz Methods
Given $A$ and $bb$ (but not $beps$), we seek to learn $bx^*$ using Kaczmarz methods. That is, we begin with an initial guess $bx^((0)) in RR^d$, and we iteratively refine it by projecting the current solution $bx^((k))$ onto problem constraints (rows of the system in @eq:orig-problem) selected by some criterion. That is, if we select row $i_k$ of the system following the $k$th iteration, we update $bx^((k))$ as follows:
#num-eq[$
  bx^((k + 1)) &= bx^((k)) + (b_i_k - innerproduct(ba_i_k, bx^((k)))) / (norm(ba_i_k)^2) ba_i_k \
  &= bx^((k)) - epsilon^((k))_i_k ba_i_k, wide "where we define" #h(0.5em) beps^((k)) = vec(epsilon^((k))_1, dots.v, epsilon^((k))_n) := A bx^((k)) - bb.
$ <eq-kaczmarz>]

#private[
    Holds because $- epsilon^((k))_i ba_i
        &= (b_i - ba_i^trp bx^((k))) ba_i$.

    Note that $epsilon_i^((k)) = ba_i^trp bx^((k)) - b_i$.

]

We hope that $bx^((k)) -> bx^*$. 
Many methods for choosing $i_k$ at each iteration have been proposed (#to-do[write about this]).

== Deriving an Optimal Selection Criterion
#let be = $bold(e)$
We define $be^((k))$ to represent the difference between the $k$th iterate and the ground truth vector $bx^*$, and $e^((k))$ to represent its squared 2-norm:
$
  be^((k)) = bx^((k)) - bx^*, wide e^((k)) = norm(be^((k)))_2^2.
$

In order to accelerate convergence of $bx^((k))$ to $bx^*$, we seek an $i_k$ at each iteration that will minimize the distance of the next estimate to the true solution, $e^((k))$. #to-do[We do this greedily, since we don't think about the whole sequence, we just think about the current iteration.]

Using the abbreviations 
#num-eq[$
    i := i_k, quad epsilon_i := epsilon_i^((k)), quad "and" quad bx := bx^((k)).
$ <eq-kth-abbreviations>]
we will expand $be^((k+1)) = norm(bx^((k+1)) - bx^*)$ using @eq-kaczmarz.

First note that
#num-eq[$
  ba_i^trp (bx - bx^*) &= (ba_i^trp bx - b_i) - (ba_i^trp bx^* - b_i) \
  &= epsilon_i - epsilon_i^*.
$ <eq-expansion-eps>]

We expand the expression for the convergence of the $(k+1)$th iterate $bx^((k+1))$ in terms of the $k$th iteration (using the abbreviations given in @eq-kth-abbreviations):
#num-eq[$
  e^((k+1)) = norm(bx^((k+1)) - bx^*)^2
  &= norm(bx - epsilon_i ba_i - bx^*)^2 \
  &= norm((bx - bx^*) - epsilon_i ba_i)^2 \
  &= norm(bx - bx^*)^2 + epsilon_i^2 norm(ba_i)^2 - 2 epsilon_i ba_i^trp (bx - bx^*) \
  &= e^((k)) + epsilon_i^2 - 2 epsilon_i (epsilon_i - epsilon_i^*) wide #[by @eq-expansion-eps] \
  &= e^((k)) - epsilon_i^2 + 2 epsilon_i epsilon_i^* \
  #private[$ 
    &= e^((k)) - epsilon_i^2 + 2 epsilon_i epsilon_i^* - (epsilon_i^*)^2 + (epsilon_i^*)^2 \
    &= e^((k)) - (epsilon_i - epsilon_i^*)^2 + (epsilon_i^*)^2 \
  $]
  &= e^((k)) - M_i
$ <eq-opt-derivation>]
where $M_i$ abbreviates $M_i^((k))$, the _optimal selection criterion_, which is defined as
#num-eq[$
  M_i^((k)) := epsilon_i^2 - 2 epsilon_i epsilon_i^* = (epsilon_i - epsilon_i^*)^2 - (epsilon_i^*)^2.
$ <eq-opt-criterion>]

Thus, if we seek to minimize the distance between the next estimate $bx$ of $bx^*$, the optimal formula for selecting $i$ at some fixed iteration is given by
#num-eq[$
  i^* = argmax_(i in {1, ..., n}) M_i, wide M_i = (epsilon_i - epsilon_i^*)^2 - (epsilon_i^*)^2.
$ <eq-i-star>]

#to-do[Write about how @eq-opt-derivation motivates maximizing $M_i$.]

== Proposing a Teacher-Optimal Selection Criterion
However, in practice, we do not have access to the true solution $bx^*$. Say we have some fixed estimate $hat(bx)$ of the true solution $bx^*$, and define its associated residual vector as 
#num-eq[$
  hat(beps) := A hat(bx) - bb.
$ <eq-teacher-residual>]
Because we hope that $hat(bx)$ can "teach" us something about $bx^*$, we call it the _teacher_. We define the _teacher-optimal selection criterion_ much as we defined the optimal selection criterion $M_i^((k))$, but using the teacher's residuals $hat(beps)$:
#num-eq[$
  hat(M)_i^((k)) :=  epsilon_i^2 - 2 epsilon_i hat(epsilon)_i = (epsilon_i - hat(epsilon)_i)^2 - hat(epsilon)_i^2,
$ <eq-teacher-opt-criterion>]
which we abbreviate as $hat(M)_i$ when the iteration index $k$ is clear or arbitrary.

To apply this criterion, at each iteration we select
#let ihat = $hat(dotless.i)$
#num-eq[$
  ihat = argmax_(i in {1, ..., n}) hat(M)_i.
$ <eq-ihat>]

In @sec-teacher-opt, we demonstrate that the teacher-optimal selection criterion $hat(M)$ (@eq-teacher-opt-criterion) converges #to-do[nicely] by showing that the optimal criterion at the teacher-optimal index, $M_(ihat)$, is not far off from $M_(i^*)$ when the teacher is reasonably accurate #to-do[(ensure this is again addressed)].

== Assumptions on Teacher
#let teachererr = $tau$
Of course, a teacher model that is wildly inaccurate will not yield near-optimal sample selection in @eq-ihat.
As we proceed in demonstrating convergence results for teacher-model online batch selection methods, we let $teachererr$ represent how near the teacher $hat(bx)$ is to the ground truth solution $bx^*$:
#num-eq[$
  teachererr := norm(hat(beps) - beps^*)_oo.
$ <eq-teacher-err>]

Thus when $teachererr = 0$, the teacher is exactly $bx^*$.

#show heading.where(level: 1): it => {
  pagebreak(weak: true)
  it
}
= Full-dataset Methods

In this section, we give convergence bounds and computational results for full-dataset teacher-based sample selection methods. That is, at each training iteration, the criterion considers _every_ row of the system as a candidate for selection, as in @eq-i-star and @eq-ihat. In practice, this is often too expensive, and we instead uniformly sample batches of rows from which to select an optimal candidate. We analyze convergence in expectation of batch-sampled methods in @sec-batch-sampled.

== Teacher-Optimal Method <sec-teacher-opt>

#to-do[
  - Theorems for teacher-optimal method (single sample)
    - Main lemmas
    - Diff eq / iterative / new method? bound proof using said lemmas
]

We demonstrate the exact relationship between the optimal criterion $M_i$ and the teacher-optimal criterion $hat(M)_i$.

#let ek = $e^((k))$
Recalling @eq-ihat, we select $ihat = argmax_i hat(M)_i$ at each iteration, so @eq-opt-derivation gives $e^((k+1)) = ek - M_ihat$. In order to bound $e^((k+1))$ from above, we then need a lower bound on $M_ihat$. We know little about the index $ihat$---it could be a clean index ($ihat in cal(C)$), or not. 

First, we pin down the exact relationship between $hat(M)_i$ and $M_i$ at an arbitrary index $i in {1, ..., n}$, using @eq-opt-criterion and @eq-teacher-opt-criterion.
#num-eq[$
  delta_i &:= hat(M)_i - M_i &= (epsilon_i^2 - 2 epsilon_i hat(epsilon)_i) - (epsilon_i^2 - 2 epsilon_i epsilon_i^*) = 2 epsilon_i (epsilon_i^* - hat(epsilon)_i).
$ <eq-delta-bound>]

Using this relationship, we may place a lower bound on $M_ihat$. Letting $j in {1, ..., n}$ represent any index, we have $hat(M)_ihat >= hat(M)_j$ by the definition of $ihat$, so
#num-eq[$
  M_ihat &= hat(M)_ihat - delta_ihat >= hat(M)_j - delta_ihat = M_j + delta_j - delta_ihat.
$ <eq-opt-at-teacher-opt>]

This bound applies for any index $j$, so we may average it for all $j$ in the clean index set $cal(C)$ to derive a bound on $M_ihat$ in terms of the clean rows of the system and a bias term.

#num-eq[$
  M_ihat >= 1/cardC sum_(j in cal(C)) M_j + 1/cardC sum_(j in cal(C)) delta_j - delta_ihat
  >= 
  underbrace(1/cardC sum_(j in cal(C)) epsilon_j^2, "perfect-teacher"\ "gain")
  - underbrace(abs(1/cardC sum_(j in cal(C))delta_j), "teacher"\ "bias on" cal(C))
  - underbrace(abs(delta_ihat), "teacher"\ "bias on" ihat).
$ <eq-bound-opt-at-teacher-opt>]

Viewing this bound in the context of @eq-opt-derivation, we have
$
  e^((k+1)) &= e^((k)) - M_ihat <= e^((k)) - 1/cardC sum_(j in cal(C)) epsilon_j^2 + abs(1/cardC sum_(j in cal(C))delta_j) + abs(delta_ihat).
$

We will further extend this bound so that it may be written in terms of only $ek$ and system constants computable before the start of training. This will allow us to apply it recursively to derive a bound on $ek$ in terms of only system constants and $e^((0))$, the initial error.

// We now further extend this bound so that it may be written in terms of system constants computable before the start of training 

=== Bounding Gain
We begin by bounding the first bracketed term in @eq-bound-opt-at-teacher-opt from below. This term is independent of the teacher, so its bound will depend only on the current iterate $bx^((k))$ (again abbreviated to $bx$) and the clean rows of the system.
#num-eq[$
  1/cardC sum_(j in cal(C)) epsilon_j^2 
  &= 1/cardC sum_(j in cal(C)) (epsilon_j - epsilon_j^*)^2 \
  &= 1/cardC norm(A_cal(C) (bx - bx^*))^2 \
  &>= (sigma_"min" (A_cal(C))^2) / cardC norm(bx - bx^*)^2 \
  &= (sigma_"min" (A_cal(C))^2) / cardC e^((k))
$ <eq-gain-bound>]

=== Bounding Teacher Bias on Clean Data

#to-do[Rewrite this section. It's a disaster.]

To bound the second term in @eq-bound-opt-at-teacher-opt, we show that #to-do[(natural language description of techniques)].

First, observe that $delta_j = -2 ba^trp_j (beps - beps^*) ba^trp_j (hat(beps) - beps^*)$ since
$
  delta_j &= 2 epsilon_j (epsilon_j^* - hat(epsilon)_j) \
  &= - 2 (epsilon_j - epsilon_j^*) (hat(epsilon)_j - epsilon_j^*) \
  &= -2 ba^trp_j (bx - bx^*) ba^trp_j (hat(bx) - bx^*).
$

This allows us to rewrite the teacher bias term as the inner product of two vectors, which we can then bound with Cauchy--Schwarz.
$
  abs(1/cardC sum_(j in cal(C))delta_j) 
  &= (2) / cardC abs(sum_(j in cal(C)) ba^trp_j (bx - bx^*) ba^trp_j (hat(bx) - bx^*)) \
  &= private(2/cardC abs((bx - bx^*)^trp (sum_(j in cal(C))  ba_j ba_j^trp) (hat(bx) - bx^*)) \ )
  &=  #purple-star\
  &= (2)/cardC abs((A_cal(C) (bx - bx^*))^trp (A_cal(C) (hat(bx) - bx^*))) \
  &<= 2/cardC norm(A_cal(C) (bx - bx^*)) norm(A_cal(C) (hat(bx) - bx^*)) wide &&#[by Cauchy--Schwarz] \
  &<= 2/sqrt(cardC) norm(A_cal(C) (bx - bx^*))_2 thin norm(A_cal(C) (hat(bx) - bx^*))_oo \
  &<= 2/sqrt(cardC) norm(A_cal(C))_2 thin norm(bx - bx^*)_2 thin norm(hat(beps) - beps^*)_oo \ 
  &<= (2 teachererr sigma_"max" (A_cal(C))) / sqrt(cardC)  sqrt(e^((k))) wide &&#[by @eq-teacher-err] \
  &<= 2 tau sqrt(e^((k))) wide &&#[because $sigma_"max" (A_cal(C)) <= norm(A_cal(C))_F = sqrt(cardC)$].
$

=== Bounding Teacher Bias on Selected Row
In bounding the last term of @eq-bound-opt-at-teacher-opt, we must first note that the teacher-optimal index $ihat$ may not lie in $cal(C)$, so we cannot make assumptions about the system noise $epsilon^*_ihat$. We instead apply the bound given in @eq-delta-bound:
$
  abs(delta_ihat) <= 2 (sqrt(e^((k))) + norm(beps^*)_oo) tau.
$

=== Combining Bounds
Putting these three bounds together into @eq-bound-opt-at-teacher-opt, we derive the desired lower bound on $M_ihat$.
$
  M_ihat &>= 
  (sigma_"min" (A_cal(C))^2) / cardC e^((k)) - 
  2 tau sqrt(e^((k))) - 
  2 (sqrt(e^((k))) + norm(beps^*)_oo) tau \
  &= (sigma_"min" (A_cal(C))^2) / cardC e^((k)) - 
  4 tau sqrt(e^((k))) - 
  2 norm(beps^*)_oo tau.
$
The error $e^((k+1))$ may thus be bounded above by a quadratic in $sqrt(ek)$:
$
  e^((k+1)) &= e^((k)) - M_ihat
  &<= (1 - sigma_"min" (A_cal(C))^2/cardC) ek + 4 tau sqrt(ek) + 2 norm(beps^*)_oo tau.
$
We state this result as @prop-teacher-opt-main-bound.

#proposition(label: <prop-teacher-opt-main-bound>)[
  When we select the training index $i_k$ using the teacher-optimal criterion (@eq-ihat), the next iteration's error $e^((k+1))$ is related to the previous error by the following inequality:
  #num-eq[$
    e^((k+1)) <= (1 - alpha) e^((k)) + beta sqrt(e^((k))) + gamma,
    \
    alpha := (sigma_"min" (A_cal(C))^2) / cardC, 
    wide beta := 4 tau, wide gamma = 2 norm(beps^*)_oo tau.
  $ <eq-master-ineq>]
]

#to-do[Explain the value of this proposition. Note that when $tau = 0$, or $hat(bx) = bx^*$, the error converges to 0 geometrically in spite of the fact that $beps^* != bold(0)$.]


=== Convergence Results
There are many ways to recursively apply @prop-teacher-opt-main-bound to bound the error at an arbitrary iteration using only constants and the initial state of training. In this section, we give three different methods, each of which gives a different bound, and none of which is strictly better than another #to-do[(justify this claim)].

==== Convergence of Teacher-Optimal Selection Criterion (Via Fixed-Point Shift)
First, we determine the fixed point of the bound in @prop-teacher-opt-main-bound. Let $r$ be the square root of this error fixed point, so that the fixed point equation becomes quadratic:
$
  r^2 = (1 - alpha) r^2 + beta r + gamma wide ==> wide r = (beta plus sqrt(beta^2 + 4 alpha gamma)) / (2 alpha).
$

(We use the upper root because $beta - sqrt(beta^2 + 4 alpha gamma) < beta - sqrt(beta^2) = 0$, but $ek$ is nonnegative.)

When $sqrt(e^((k))) = r$, @prop-teacher-opt-main-bound gives $sqrt(e^((k+1))) <= sqrt(e^((k)))$; that is, the bound does not loosen or tighten.

#let dk = $d^((k))$
Assume $sqrt(e^((k))) >= r$, and define $dk := sqrt(e^((k))) - r$ to be the distance of the square root error from the fixed point, and rewrite @prop-teacher-opt-main-bound as a quadratic in $d^((k))$ using $e^((k)) = (dk + r)^2$. 
As we have done with $bx^((k))$, we abbreviate $d^((k))$ to $d$.
#num-eq[$
  (d^((k+1)) + r)^2 = e^((k+1)) &<= (1 - alpha) (d + r)^2 + beta (d + r) + gamma \
  &= (1 - alpha) (d^2 + 2  d r + r^2) + beta (d + r) + gamma \
  &= (1 - alpha) d^2 + (2 r (1 - alpha) + beta) d + underbrace((1 - alpha) r^2 + beta r + gamma, r^2) private(\
  &<=(rho d + r)^2\, wide rho := max(sqrt(1 - alpha), thick 1 - alpha + beta/(2 r)).)
$ <eq-fixed-shift-bound>]

Expanding both sides of @eq-fixed-shift-bound, we may derive a bound on $d^((k+1))$ in terms of $d$, as follows:
#num-eq[$
(d^((k+1)) + r)^2 &<= (1 - alpha) d^2 + (2 r (1 - alpha) + beta) d + r^2 \
==> d^((k+1)) &<= sqrt((1 - alpha) d^2 + (2 r (1 - alpha) + beta) d + r^2) - r.
$]

We seek a constant $rho$ such that $d^((k+1)) <= rho d$. By assuming $d > 0$, the above inequality imposes a satisfactory $rho$.
$
  d^((k+1)) = sqrt((1 - alpha) d^2 + (2 r (1 - alpha) + beta) d + r^2) - r &<= rho d \
  <==> (1 - alpha) d^2 + (2 r (1 - alpha) + beta) d + r^2 &<= rho^2 d^2 + 2rho r d + r^2
$
$
  <==> 1 - alpha <= rho^2 quad "and" quad 2 r(1-alpha) + beta <= 2 rho r \
  <==> rho >= sqrt(1 - alpha) quad "and" quad rho >= 1 - alpha + beta / (2 r) \
  <== rho := max(sqrt(1 - alpha), thick thick 1 - alpha + beta / (2 r))
$

Recursively applying $d^((k+1)) <= rho d$, we may bound $d$ in terms of $d^((0))$:
#num-eq[$ d <= rho d^((k-1)) <= dots.c <= rho^k d^((0)). $ <eq-d-bound-recursion>]

Thus $dk = sqrt(e^((k))) - r$ converges geometrically toward 0. Rewriting this in terms of $ek$, we may add $r$ to both sides of @eq-d-bound-recursion and square to derive
$
  (dk + r)^2 &<= (rho^k d^((0)) + r)^2,
$
which, by application of the definition of $dk$, gives
$
  ek &<= (rho^k sqrt(e^((0)) - r) + r)^2.
$

==== Convergence of Teacher-Optimal Selection Criterion (Via Peter--Paul)
@prop-teacher-opt-main-bound gives a bound on $e^((k+1))$ in terms of $e^((k))$. We now use it to derive a bound on $e^((k))$ in terms of only $e^((0))$ and other information available from the start of training.

#let efix = $e_"fixed"$
#theorem(label: <thm-master>)[
  Using the teacher-optimal criterion (@eq-ihat), the squared #to-do[(make this more precise)] error converges geometrically toward 
  $
    efix = 2 (beta^2/alpha^2 + gamma / alpha),
  $
  with
  $
    e^((k)) <= a^k (e^((0)) - efix) + efix.
  $ 
  where $a := 1 - alpha\/2$.
  #to-do[reorder things, bound #efix]
]


==== Convergence of Teacher-Optimal Selection Criterion (Via Differential Equation)
To estimate $partial / (partial k) e^((k))$, we can use the forward Euler method with step size 1:
$
  partial/(partial k) e^((k)) approx (e^((k+1)) - e^((k))) / 1.
$
#to-do[Explain why it's actually not just $approx$. It's actually an inequality once $k$ is large enough]

By @eq-master-ineq, we may write

$
  e^((k+1)) &= (1 - alpha) e^((k)) + beta sqrt(e^((k))) + gamma \
  ==> e^((k+1)) - e^((k)) &<= f(e^((k)))
$ 
where
$
  f(ek) := -alpha ek + beta sqrt(ek) + gamma,
$

so we may write the bound
#let dif = $d$
$
  dif/(dif k) ek <= f(ek)
$

Define $dif y \/ dif t = f(y) = -alpha y + beta sqrt(y) + gamma$, and note that when $k = t$,
$
  dif / (dif k) e^((k)) <= (dif y) / (dif t) = f(y),
$

so $e^((k)) <= y(t)$ when $k = t$. 

We begin by seeking an error floor given by this bound. This will occur when $dif y \/ dif t = 0$.
$
  0 &= (dif y) / (dif t) \
  ==> 0 &= -alpha y + beta sqrt(y) + gamma \
  ==> 0 &= -alpha z^2 + beta z + gamma wide #[where $z = sqrt(y)$] \
  ==> z &= (-beta - sqrt(beta^2 + 4 alpha gamma)) / (-2 alpha) \
  ==> y &= ((-beta - sqrt(beta^2 + 4 alpha gamma)) / (-2 alpha))^2 =: y_"floor"
$

Thus, given some error $macron(e) > y_"floor"$, if $y(t) <= macron(e)$, certainly $e^((k)) <= macron(e)$ as well.

Fix $macron(e) > y_"floor"$. We seek the $T$ such that $y(t) = macron(e)$, thereby giving an upper bound on the lowest $k$ such that $e^((k)) <= macron(e)$.

$
  (dif y) / (dif t) &= f(y(t)) \
  ==> dif t &= (dif y) / (f(y(t))) \
  ==> T &= integral_(y(0))^(y(T)) (dif y) / (f(y(t))) \
  &= integral_(e^((0)))^macron(e) (dif y) / (-alpha y + beta sqrt(y) + gamma) \
  &= (2 / D) (r_+ log (sqrt(e^((0))) - r_+) / (sqrt(macron(e)) - r_+) - r_- log (sqrt(e^((0))) - r_-) / (sqrt(macron(e)) - r_-)),
$
$
  D := sqrt(beta^2  + 4 alpha gamma), wide r_plus.minus = (beta plus.minus D) / (2 alpha)
$

#to-do[Refactor into theorem format.]

#to-do[Explain why the integral holds. #link("https://chatgpt.com/share/e/6ace6d75-7e20-83e8-a20a-af4b368e6f7d")] 

== Reducible Holdout Loss Selection
In #link("https://arxiv.org/pdf/2206.07137")[_Prioritized Training on Points that are Learnable, Worth Learning, and Not Yet Learnt_], Mindermann et al. propose (in the context of deep learning) another sample selection method, _reducible holdout loss selection (RHO-LOSS)._ We may use this method in a modified form #footnote[Mindermann et al. define RHO-LOSS using a "holdout" dataset, which is simply a subset of the training data reserved for training the teacher model (here, our teacher model is $hat(bx)$). Because we assume that our teacher is an arbitrary fixed vector and do not discuss the "training" of our teacher, we instead define the holdout loss using the teacher residual.] for the Kaczmarz method. Using the same notation as in @sec-introduction, we define RHO-LOSS by the selection criterion
$
  hat(R)_i := epsilon_i^2 - hat(epsilon)_i^2, wide ihat = argmax_(i in {1, ..., n}) hat(R)_i.
$

We begin with an approach similar to that used for the teacher-optimal method in @sec-teacher-opt.
Again, @eq-opt-derivation gives $e^((k+1)) = ek - M_ihat$, so we seek a lower bound on $M_ihat$ in order to prove convergence of $ek$ to within some neighborhood of $0$.
We begin by characterizing the relationship between $hat(R)_i$ and $M_i$ at an arbitrary index $i in {1, ..., n}$ (compare @eq-delta-bound):

#num-eq[$
  delta_i := hat(R)_i - M_i = (epsilon_i^2 - hat(epsilon)_i^2) - (epsilon_i^2 - 2 epsilon_i epsilon_i^*) = 2 epsilon_i epsilon_i^* - hat(epsilon)_i^2. 
$ <eq-delta-bound-rho>]

Using this relationship to begin to bound $M_ihat$, we find that for any index $j in {1, ..., n}$ (compare @eq-opt-at-teacher-opt),
#num-eq[$
  M_ihat = hat(R)_ihat - delta_ihat >= hat(R)_j - delta_ihat = M_j + delta_j - delta_ihat.
$ <eq-opt-at-rho-opt>]

Because @eq-opt-at-rho-opt holds for any $j$, we may average each side of the inequality over every $j in cal(C)$ as before (compare @eq-bound-opt-at-teacher-opt):
#num-eq[$
  M_ihat >= 1/cardC sum_(j in cal(C)) M_j + 1/cardC sum_(j in cal(C)) delta_j - delta_ihat
  >= 
  underbrace(1/cardC sum_(j in cal(C)) epsilon_j^2, "perfect-teacher"\ "gain")
  - underbrace(abs(1/cardC sum_(j in cal(C))delta_j), "teacher"\ "bias on" cal(C))
  - underbrace(abs(delta_ihat), "teacher"\ "bias on" ihat).
$ <eq-bound-opt-at-rho>]

Note that this inequality is identical to that in @eq-bound-opt-at-teacher-opt, except in that $delta_i$ is now defined by @eq-delta-bound-rho.

=== Bounding Gain
The lower bound on the first bracketed term is given by @eq-gain-bound.

=== Bounding Teacher Bias on Clean Data
We now give an upper bound on the second term in @eq-bound-opt-at-rho. #to-do[Complete after rewriting analogue in teacher-opt section.]

// === Bounding Teacher Bias on Selected Row
// In bounding the last term of @eq-bound-opt-at-teacher-opt, we must first note that the teacher-optimal index $ihat$ may not lie in $cal(C)$, so we cannot make assumptions about the system noise $epsilon^*_ihat$. We instead apply the bound given in @eq-delta-bound:
// $
//   abs(delta_ihat) <= 2 (sqrt(e^((k))) + norm(beps^*)_oo) tau.
// $




#to-do[
  #v(10em)
  - Theorems for RHO-LOSS (single sample)
    - Modification of teacher-optimal lemmas for RHO-LOSS if possible?
    - Bound proof
  ]

== Empirical Results
#to-do[
  - Computational results for single-sample methods, comparison to (probably very loose) bound
    - Gaussian system
    - Real datasets
    - Clock-time comparison to Randomized Kaczmarz? RK should win on large datasets probably? Then batch-sampled opt method should blow RK out of the water
]
= Batch-sampled Methods <sec-batch-sampled>

teststets
#to-do[
  - Proof of convergence in expectation for teacher-optimal method
    - Lemmas
    - Proof from lemmas (unless the proof is too complex for that?)

  - Proof of convergence in expectation for RHO-LOSS also?

  - Computational results for batch-sampled teacher-based Kaczmarz
]

= Appendix

#let proof(..args) = {
  theoretic.presets.basic.proof(suffix: sym.square, ..args)
  divider()
}

#to-do[
  - Prove lemmas for each section
]

#proof(title: [@thm-master])[
  We begin with an application of the Peter-Paul inequality to absorb the $sqrt(e^((k)))$ term of @eq-master-ineq into the $e^((k))$ term. For any $eta > 0$,
  $
    beta sqrt(e^((k)))  &<= beta^2/(2 eta) + (eta e^((k)))/2.
  $
  To ensure that $1 - alpha + eta\/2 < 1$, we set $eta = alpha$, which gives the following relaxation of @eq-master-ineq in terms of only $e^((k))$ and constants:
  #num-eq[$
    e^((k+1)) &<= private((1 - alpha/2) e^((k)) + beta^2 / alpha + gamma &=) a e^((k)) + c, wide a := 1 - alpha/2, quad c := beta^2/alpha + gamma.
  $ <eq-relaxed-master>]
  
  Note that the fixed point #efix of this inequality is 
  $
    efix = 2 ((beta^2)/(alpha^2) + gamma/alpha),
  $
  since
  $
    efix private(= efix - (beta^2 / alpha + gamma) + beta^2 / alpha + gamma) = a efix + c,
  $
  and it follows that 
  $
    efix = c / (1 - a).
  $
  
  We now derive a bound on $e^((k))$ in terms of only $e^((0))$ and constants. By repeated application of @eq-relaxed-master, we have
  $
    e^((k)) <= a e^((k-1)) + c <= a (a e^((k-2)) + c) + c <= dots.c &<= a^k e^((0)) + sum_(ell=0)^(k-1) a^ell c.
  $
  By the finite geometric series formula, $sum_(ell=0)^(k-1) a^ell c = (c (1 - a^k)) / (1 - a)$, so
  $
    e^((k)) &<= a^k e^((0)) + (c (1 - a^k)) / (1 - a) \
    private(&= a^k e^((0)) + c/(1 - a) - (c) / (1 - a) a^((k)) \ )
    &= a^k (e^((0)) - efix) + efix.
  $
]