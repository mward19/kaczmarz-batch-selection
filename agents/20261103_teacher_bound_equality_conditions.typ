#import "@local/quicknotes:0.1.1": *
#show: template.with(font-mode: "computer")

#set page(paper: "a4", margin: 22mm, numbering: "1")
#set text(font: "New Computer Modern", size: 10.5pt)
#set par(justify: true, leading: 0.6em)
#set heading(numbering: "1.1")
#set math.equation(numbering: "(1)")
#let trp = $"T"$
#let bx = $bold(x)$
#let ba = $bold(a)$
#let bb = $bold(b)$
#let beps = $bold(epsilon)$
#let be = $bold(e)$
#let cardC = $abs(cal(C))$
#let ihat = $hat(dotless.i)$

#align(center)[
  #text(size: 17pt, weight: "bold")[Equality conditions for the teacher-optimal bound] \
  GPT-6 (Codex) \
  November 3, 2026
]

= Scope and conclusion

This report audits the inequalities used to obtain `prop-teacher-opt-main-bound` in `paper.typ`, including the supporting lemmas and their appendix proofs. Later convergence relaxations are outside this proposition's proof.

When $e^((k))>0$ and $tau>0$, no system satisfying the assumptions below makes the entire proof tight. The signed equality conditions for the clean-row and selected-row teacher biases conflict. In addition, the final spectral relaxation would require $A_cal(C)$ to have rank one. Exact equality is possible when the teacher is perfect or the current iterate is already correct; examples follow.

= Notation and assumptions

Fix an iteration $k$, and use the paper's notation. The iterate error vector is $be^((k))=bx^((k))-bx^*$, with $e^((k))=norm(be^((k)))_2^2$. The teacher error vector is $hat(bx)-bx^*$, and $tau=norm(hat(beps)-beps^*)_oo$ as in `eq-teacher-err`. For every row $i$,
$
 epsilon_i^((k))-epsilon_i^* = ba_i^trp be^((k)), quad
 hat(epsilon)_i-epsilon_i^* = ba_i^trp (hat(bx)-bx^*), \
 M_i^((k))=(epsilon_i^((k))-epsilon_i^*)^2-(epsilon_i^*)^2, \
 delta_i=2epsilon_i^((k))(epsilon_i^* - hat(epsilon)_i), quad
 hat(M)_i^((k))=M_i^((k))+delta_i.
$
The selected index is $ihat$, as in `eq-ihat`. The exact Kaczmarz update gives
$
 e^((k+1))=e^((k))-M_ihat^((k)).
$
We assume real matrices, unit rows, $n>d$, $cardC>d$, $epsilon_j^*=0$ for $j in cal(C)$, and full column rank of $A_cal(C)$. The paper states that $cardC>d$ implies full column rank, but that implication is false; full rank is a separate required assumption. It makes $tau=0$ equivalent to $hat(bx)=bx^*$.

All equality conditions below concern a single iteration. Equality for one step need not persist along a trajectory.

= Audit of the selected-row residual bound

The appendix proof of `eq-delta-bound` bounds $abs(delta_i)$. Each step can be tight individually, with the following conditions.

+ *Triangle inequality:* 
  $abs(epsilon_i^((k))) <= abs(epsilon_i^((k))-epsilon_i^*)+abs(epsilon_i^*)$.
  Equality holds exactly when $(epsilon_i^((k))-epsilon_i^*)epsilon_i^*>=0$.

+ *Row Cauchy–Schwarz:* 
  $abs(ba_i^trp be^((k)))<=norm(be^((k)))_2=sqrt(e^((k)))$.
  If $e^((k))>0$, equality holds exactly when $ba_i=plus.minus be^((k))/sqrt(e^((k)))$; if $e^((k))=0$, it is automatic.

+ *Noise maximum:* $abs(epsilon_i^*)<=norm(beps^*)_oo$. Equality means $abs(epsilon_i^*)=norm(beps^*)_oo$.

+ *Absolute product:* $abs(delta_i)=2abs(epsilon_i^((k)))abs(hat(epsilon)_i-epsilon_i^*)$ is an identity.

+ *Teacher maximum:* $abs(hat(epsilon)_i-epsilon_i^*)<=tau$. In the product above, this step is tight if this absolute value is $tau$ or $epsilon_i^((k))=0$.

+ *Combined estimate:* if $tau>0$ and $sqrt(e^((k)))+norm(beps^*)_oo>0$, equality in
  $abs(delta_i)<=2(sqrt(e^((k)))+norm(beps^*)_oo)tau$
  holds exactly when
  $
    abs(ba_i^trp be^((k)))=sqrt(e^((k))), quad
    abs(epsilon_i^*)=norm(beps^*)_oo, \
    (ba_i^trp be^((k)))epsilon_i^*>=0, quad
    abs(hat(epsilon)_i-epsilon_i^*)=tau.
  $
  If either multiplying factor on the right is zero, both sides vanish. A one-dimensional unit row can attain this residual estimate, but need not be selected by the teacher criterion.

= Selection, averaging, and signs

+ *Teacher maximization:* for any comparison row $j$,
  $
    M_ihat^((k))=hat(M)_ihat^((k))-delta_ihat
    >=hat(M)_j^((k))-delta_ihat
    =M_j^((k))+delta_j-delta_ihat.
  $
  Equality holds exactly when $hat(M)_j^((k))=hat(M)_ihat^((k))$.

+ *Averaging over clean rows:* since $epsilon_j^*=0$ on $cal(C)$,
  $
    M_ihat^((k))>=
    1/cardC sum_(j in cal(C)) (ba_j^trp be^((k)))^2
    + 1/cardC sum_(j in cal(C)) delta_j-delta_ihat.
  $ <signed-average>
  Equality holds exactly when every clean row ties the selected row in teacher score: $forall j in cal(C), hat(M)_j^((k))=hat(M)_ihat^((k))$. The gap is the average of the nonnegative differences between those teacher scores.

+ *Clean-bias absolute value:* replacing the clean average bias by its negative absolute value is tight exactly when $sum_(j in cal(C))delta_j<=0$.

+ *Selected-bias absolute value:* replacing $-delta_ihat$ by $-abs(delta_ihat)$ is tight exactly when $delta_ihat>=0$, or equivalently
  $(epsilon_ihat^((k)))(hat(epsilon)_ihat-epsilon_ihat^*)<=0$.

These sign conditions are essential to simultaneous equality.

= Gain and clean teacher bias

+ *Smallest singular value:* the gain bound is
  $
    1/cardC sum_(j in cal(C)) (ba_j^trp be^((k)))^2
    =norm(A_cal(C)be^((k)))_2^2/cardC
    >=alpha e^((k)),
  $
  where $alpha=sigma_"min"(A_cal(C))^2/cardC$ is defined in `eq-master-ineq`. If $e^((k))>0$, equality holds exactly when $be^((k))$ is a right singular vector for the smallest singular value. If $e^((k))=0$, equality is automatic.

+ *Clean Cauchy–Schwarz:* 
  $
  abs((A_cal(C)be^((k)))^trp A_cal(C)(hat(bx)-bx^*))
  <=norm(A_cal(C)be^((k)))_2 norm(A_cal(C)(hat(bx)-bx^*))_2.
  $
  Equality holds if either vector is zero or they are real scalar multiples. Under full column rank, if $e^((k))>0$ and $tau>0$, the teacher error must be a nonzero scalar multiple of $be^((k))$. Equality with the required nonnegative sign of the clean bias requires a positive multiple.

+ *Euclidean and maximum norms:* 
  $norm(A_cal(C)(hat(bx)-bx^*))_2<=sqrt(cardC)norm(A_cal(C)(hat(bx)-bx^*))_oo$.
  Equality holds exactly when all clean-row teacher residual errors have equal absolute value. In a product with $norm(A_cal(C)be^((k)))_2=0$, the product vanishes without requiring equality here.

+ *Largest singular value and clean restriction:* 
  $norm(A_cal(C)be^((k)))_2<=sigma_"max"(A_cal(C))sqrt(e^((k)))$
  and
  $norm(A_cal(C)(hat(bx)-bx^*))_oo<=tau$.
  For nonzero factors, equality requires respectively a largest right singular vector and a clean row attaining the global teacher-error maximum. Equality in the preceding maximum-norm comparison additionally requires every clean row to attain that maximum.

+ *Simultaneous smallest and largest singular values:* for $e^((k))>0$, simultaneous equality in the gain bound and the largest-singular-value bound requires
  $
  sigma_"min"(A_cal(C))=sigma_"max"(A_cal(C)), quad
  A_cal(C)^trp A_cal(C)= (cardC/d) I_d.
  $ <tight-frame>
  The second identity follows from full column rank and row normalization, which gives $op("tr")(A_cal(C)^trp A_cal(C))=cardC$.

+ *Final spectral relaxation:* 
  $sigma_"max"(A_cal(C))<=norm(A_cal(C))_F=sqrt(cardC)$.
  Equality in the inequality holds exactly when $A_cal(C)$ has rank one. Under full column rank this means $d=1$. In the bias product the relaxation is also tight when $tau sqrt(e^((k)))=0$, since then the product vanishes regardless of the matrix inequality.

The primitive clean norm bounds can all be attained at once: take $d=1$, every clean row equal to $1$, and positive scalar iterate and teacher errors. That example does not satisfy the selected-bias sign condition at a row saturating `eq-delta-bound`. For $d>1$, repeated coordinate rows attain the smallest and largest singular-value bounds together, but the Frobenius relaxation is strict when $tau>0$ and $e^((k))>0$.

= Simultaneous equality when both errors are positive

Suppose $e^((k))>0$ and $tau>0$. Equality in the proposition's derivation requires clean rows to tie the selected teacher score, the clean average bias to be nonpositive, and $delta_ihat>=0$. It also requires simultaneous equality in the gain and largest-singular-value bounds, equality in clean Cauchy–Schwarz, all clean teacher residual errors to attain absolute value $tau$, and equality in `eq-delta-bound` at $ihat$.

Already the signed conditions make this impossible. Equality in clean Cauchy–Schwarz with the clean bias sign requires
$hat(bx)-bx^*=lambda be^((k))$
for $lambda>0$. Equality in `eq-delta-bound` at $ihat$ gives
$abs(ba_ihat^trp be^((k)))=sqrt(e^((k)))$
and
$(ba_ihat^trp be^((k)))epsilon_ihat^*>=0$.
Therefore
$
 (epsilon_ihat^((k)))(hat(epsilon)_ihat-epsilon_ihat^*)
 = (ba_ihat^trp be^((k)) + epsilon_ihat^*)(lambda ba_ihat^trp be^((k)))>0,
$
so $delta_ihat<0$. This contradicts the selected-bias equality requirement $delta_ihat>=0$. Since each relaxation contributes a nonnegative gap to the lower bound for $M_ihat^((k))$, the contradiction proves that the final inequality `eq-master-ineq` is strict whenever $e^((k))>0$ and $tau>0$.

== A near-sharp family with positive iterate and teacher errors

Although equality is impossible in this regime, the bound can be arbitrarily close in relative terms. Let $R>0$ and $s>0$, and take $d=2$, $n=5$, with rows and data
 $
 A=mat(1,0; 1,0; 0,1; 0,1; 1,0), quad
 bx^*=vec(0,0), quad
 bb=vec(0,0,0,0,2s),
 $
where the first four rows are clean and the fifth is corrupted. In the paper's residual convention $beps^*=A bx^*-bb$, so
 $
 beps^*=vec(0,0,0,0,-2s), quad
 hat(bx)=vec(s,s), quad
 bx^((k))=vec(R,R).
 $
Thus $tau=norm(A(hat(bx)-bx^*))_oo=s>0$, $B:=norm(beps^*)_oo=2s$, and $e^((k))=2R^2>0$. The clean frame has $A_cal(C)^trp A_cal(C)=2I_2$, hence $alpha=2/4=1/2$.

The current residual and teacher residual vectors are
 $
 beps^((k))=A bx^((k))-bb=vec(R,R,R,R,R-2s), \\
 hat(beps)=A hat(bx)-bb=vec(s,s,s,s,-s).
 $
Every clean row has teacher score $hat(M)_j^((k))=R^2-2R s$. The fifth row ties it because
 $
 hat(M)_5^((k))=(R-2s)^2-2(R-2s)(-s)=R^2-2R s.
 $
Therefore choosing the fifth row is permitted by the argmax rule (the example uses a tie resolved in favor of that row). Its Kaczmarz update is $bx^((k+1))=vec(2s,R)$, and hence
 $
 e^((k+1))=R^2+4s^2.
 $
The right-hand side of `eq-master-ineq` is
 $
 (1-alpha)e^((k))+4tau sqrt(e^((k)))+2B tau
 =R^2+4sqrt(2)R s+4s^2.
 $
So the inequality is strict, with additive gap $4sqrt(2)R s$, but
 $
 frac{e^((k+1))}{(1-alpha)e^((k))+4tau sqrt(e^((k)))+2B tau}
 =frac{R^2+4s^2}{R^2+4sqrt(2)R s+4s^2} -> 1 quad
 "as" quad s/R -> 0.
 $
This is a meaningful worst-case sense in which the bound is close: for every member of the family both errors are nonzero, yet the relative slack can be made arbitrarily small. The closeness is multiplicative; the additive slack itself is not zero for finite $R,s>0$. The construction also makes the dependence on tie-breaking explicit: if a clean row is selected instead, its update gives $e^((k+1))=R^2$, below the displayed near-worst-case value.

= Degenerate cases and exact equality examples

== Perfect teacher: $tau=0$

Full column rank implies $hat(bx)=bx^*$. Equality in `eq-master-ineq` holds exactly when
$
 A_cal(C)^trp A_cal(C)be^((k))=sigma_"min"(A_cal(C))^2 be^((k)), \
 forall j in cal(C), (ba_j^trp be^((k)))^2=alpha e^((k)), \
 max_(1<=i<=n) ((ba_i^trp be^((k)))^2-(epsilon_i^*)^2)=alpha e^((k)).
$
When $e^((k))=0$, the first condition is automatic and the maximum is zero because there is a clean row. For positive error, these conditions say the iterate error is in the smallest singular subspace, every clean row attains the same gain, and no row has a larger criterion.

An exact example in any dimension uses two copies of each coordinate row, $n=cardC=2d$, all rows clean, $bx^*=hat(bx)=bb=beps^*=bold(0)$, and $bx^((k))=q(1,dots,1)^trp$ for $q>0$. Then $A_cal(C)^trp A_cal(C)=2I_d$, $alpha=1/d$, and every row criterion is $q^2$. Every selected update leaves squared error $(d-1)q^2=(1-alpha)e^((k))$, attaining the proposition.

== Correct iterate and imperfect teacher: $e^((k))=0$, $tau>0$

The proposition becomes $e^((k+1))<=2norm(beps^*)_oo tau$. If $norm(beps^*)_oo=0$, all current residuals vanish, all criteria tie at zero, and equality holds for any selection.

If $norm(beps^*)_oo>0$, equality holds exactly when the selected row satisfies
$
 abs(epsilon_ihat^*)=norm(beps^*)_oo=2tau, \
 abs(hat(epsilon)_ihat-epsilon_ihat^*)=tau, \
 epsilon_ihat^*(hat(epsilon)_ihat-epsilon_ihat^*)=-norm(beps^*)_oo tau,
$
and every row has $hat(M)_i^((k))<=0$. Then the clean rows and the selected row tie at zero. Necessity follows from equality in the selected-row residual bound and the selected-bias sign step; conversely these conditions give equality directly.

For an example, let $d=1$, $n=3$, $cal(C)={1,2}$, $A=mat(1;1;1)$, $bx^*=bx^((k))=0$, $hat(bx)=T$ for $T>0$, $beps^*=vec(0,0,-2T)$, and $bb=vec(0,0,2T)$. Then $tau=T$, $norm(beps^*)_oo=2T$, and all teacher criteria are zero. Selecting the third row is allowed; the update gives $e^((k+1))=4T^2=2norm(beps^*)_oo tau$. Thus exact equality can depend on how a tie is resolved.

These cases cover all possibilities for whether $e^((k))$ and $tau$ vanish.

= Bounds that keep the paper's notation

== Keep the clean largest singular value

Omitting the final spectral-to-Frobenius relaxation gives
$
 M_ihat^((k))>=alpha e^((k))
 -2tau(1+sigma_"max"(A_cal(C))/sqrt(cardC))sqrt(e^((k)))
 -2norm(beps^*)_oo tau,
 \
 e^((k+1))<=(1-alpha)e^((k))
 +2tau(1+sigma_"max"(A_cal(C))/sqrt(cardC))sqrt(e^((k)))
 +2norm(beps^*)_oo tau.
$ <spectral-improvement>
For full-rank $A_cal(C)$ with $d>1$, $sigma_"max"(A_cal(C))<sqrt(cardC)$, so this improves the displayed bound in `eq-master-ineq` when $tau sqrt(e^((k)))>0$. Also,
$
 sigma_"max"(A_cal(C))^2/cardC<=1-(d-1)alpha,
$
because the other $d-1$ squared singular values are at least $sigma_"min"(A_cal(C))^2$.

== Preserve the clean-row signs

The selection rule itself yields the exact lower bound
$
 M_ihat^((k)) >=
 max_(j in cal(C)) ((ba_j^trp be^((k)))^2
 -2(ba_j^trp be^((k))) (ba_j^trp(hat(bx)-bx^*))) \
 +2(epsilon_ihat^((k)))(hat(epsilon)_ihat-epsilon_ihat^*).
$ <signed-improvement>
The first comparison is tight when the selected row ties a maximizing clean row; averaging the clean comparisons recovers the proof's clean-gain term. This form keeps the selected sign information that the absolute-value estimates discard.

The exact update also gives the teacher-independent bound
$
 e^((k+1))=e^((k))-(ba_ihat^trp be^((k)))^2+(epsilon_ihat^*)^2
<=e^((k))+norm(beps^*)_oo^2.
$
One can take the minimum of this and `eq-master-ineq`.

== A sharp one-dimensional bound

In $d=1$, every normalized row is $plus.minus 1$. Thus $(ba_i^trp be^((k)))^2=e^((k))$ and
$(ba_i^trp be^((k)))(ba_i^trp(hat(bx)-bx^*))=be^((k))(hat(bx)-bx^*)$ for every row. Comparing the selected score to a clean row gives
$
 (epsilon_ihat^*)^2+2epsilon_ihat^*(hat(epsilon)_ihat-epsilon_ihat^*)<=0,
 quad abs(epsilon_ihat^*)<=2tau.
$
The one-dimensional update leaves error $(epsilon_ihat^*)^2$, so
$
 e^((k+1))<=(norm(beps^*)_oo)^2, quad e^((k+1))<=4tau^2.
$ <sharp-one-dimensional>
Hence $e^((k+1))<=min((norm(beps^*)_oo)^2,4tau^2)$, independent of the current iterate error. This is sharp: the preceding three-row example with $bx^((k))=R>0$ instead of zero has all teacher scores tied and attains $e^((k+1))=4T^2$. More generally, choosing the selected row's noise magnitude below or above $2tau$ attains respectively the first or second term of the minimum, with an additional row if needed to set $norm(beps^*)_oo$.

= Corrections to the paper's derivation

The clean-bias paragraph begins with an expression involving $ba_j^trp(epsilon-epsilon^*)$. These are incompatible dimensions: residual vectors have length $n$, while $ba_j$ has length $d$. The intended vectors are $bx-bx^*$ and $hat(bx)-bx^*$.

The replacement of $epsilon_j$ by $epsilon_j-epsilon_j^*$ is valid because $j$ is clean. State that restriction explicitly; for a general row the noise term remains. Finally, specify whether the selection rule permits arbitrary maximizers or uses a fixed tie rule. Equality for an imperfect teacher requires a corrupted row tied with every clean row, so a sharpness claim should state its dependence on tie resolution.
