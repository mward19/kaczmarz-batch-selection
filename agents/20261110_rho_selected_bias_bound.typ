#import "@local/quicknotes:0.1.1": *
#show: template.with(font-mode: "computer")

#set page(paper: "a4", margin: 22mm, numbering: "1")
#set text(font: "New Computer Modern", size: 10.5pt)
#set par(justify: true, leading: 0.6em)
#let beps = $bold(epsilon)$
#let ihat = $hat(dotless.i)$

#align(center)[
  #text(size: 17pt, weight: "bold")[A Bound on the RHO-LOSS Selected-Row Bias] \
  GPT-6.1 Sol \
  November 10, 2026
]

= Setup

Assume the rows are normalized and that
$norm(hat(beps)-beps^*)_2<=tau$.
Then $abs(epsilon_i-epsilon_i^*)<=sqrt(e^((k)))$, $abs(epsilon_i^*)<=norm(beps^*)_oo$, and $abs(hat(epsilon)_i-epsilon_i^*)<=tau$.

= Exact decomposition

Starting from $delta_i=2epsilon_i epsilon_i^* - hat(epsilon)_i^2$, expand around $epsilon_i^*$:
$
  delta_i
  &= 2(epsilon_i-epsilon_i^*)epsilon_i^*
     - (hat(epsilon)_i-epsilon_i^*)^2
     - 2epsilon_i^*(hat(epsilon)_i-epsilon_i^*)
     + (epsilon_i^*)^2.
$
The cross term has a minus sign, but its value can be positive or negative because it depends on the signs of $epsilon_i^*$ and $hat(epsilon)_i-epsilon_i^*$.

= One-sided bound for the selected row

For the selected-row contribution, the comparison inequality contains $-delta_ihat$. It is enough to upper-bound $delta_i$. Since $-(hat(epsilon)_i-epsilon_i^*)^2<=0$ and
$-2epsilon_i^*(hat(epsilon)_i-epsilon_i^*)<=2abs(epsilon_i^*)abs(hat(epsilon)_i-epsilon_i^*)$,
$
  delta_i
  &<= 2abs(epsilon_i-epsilon_i^*)abs(epsilon_i^*)
     +2abs(epsilon_i^*)abs(hat(epsilon)_i-epsilon_i^*)
     +(epsilon_i^*)^2 \
  &<= 2sqrt(e^((k)))norm(beps^*)_oo
     +2tau norm(beps^*)_oo
     +norm(beps^*)_oo^2.
$
Therefore
$
  -delta_ihat
  &>= -2sqrt(e^((k)))norm(beps^*)_oo
      -2tau norm(beps^*)_oo
      -norm(beps^*)_oo^2.
$
The negative square can be discarded for this upper bound; it cannot in general be replaced by $-norm(hat(beps)-beps^*)_oo^2$, since a particular coordinate need not attain the vector's maximum norm. Likewise, the cross term cannot be bounded above by a negative norm term without an additional sign assumption.

If instead $tau$ denotes an RMS coordinate error and $norm(hat(beps)-beps^*)_2<=sqrt(n)tau$, replace each $tau$ above by $sqrt(n)tau$. With the paper's Euclidean teacher-error definition, no $sqrt(n)$ factor is needed.
