#import "oxthesis.typ": *

#show: oxthesis.with(
  title: " Predicting cancer with Bayesian monotonic spline models",
  author: "Stanley E. Lazic",
  abstract: include "abstract.typ",
  acknowledgements: include "acknowledgements.typ",
  font-size: 12pt
)

#include "introduction.typ"

#pagebreak()
#include "methodology.typ"

#pagebreak()
#include "simulations.typ"

#pagebreak()
#include "results.typ"

#pagebreak()
#include "discussion.typ"

#pagebreak()
#include "references.typ"

#counter(heading).update(0)

#pagebreak()
#include "appendix.typ"
