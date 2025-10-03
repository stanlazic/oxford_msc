#import "@preview/wordometer:0.1.4": word-count, total-words, word-count-of, string-word-count
#show: word-count


// parse-font() // {{{
#let parse-font(font-size) = {
  if font-size < 10pt {
    panic([Font font-size must be at least 10pt])
  } else {
    font-size
  }
}
//  // }}}

// init-title-page() // {{{
#let init-title-page(title,
                author
              ) = {
  set par(spacing: 0em)

  set page(
    paper: "a4"
    // width: width,    
    // height: height,
    // margin: (
    //   top: top-margin,
    //   bottom: bottom-margin,
    //   left: 32mm,
    //   right: 32mm
    // ),
  ) 

  align(center)[
    #v(12%)
    #text(size: 1.4em)[#title]
    #v(12%)
    #image("img/beltcrest.svg")
    #v(12%)
    #author
    #v(20%)
    This dissertation has been submitted to the University of Oxford in partial fulfilment of the
requirement for the award of the degree of MSc in Evidence-Based Health Care (Medical
Statistics).
    #v(10%)
	]
	align(left)[
		Date: 2026 \
		Candidate Number: XXX
	]

}
//  // }}}

// init-abstract() // {{{
#let init-abstract(abstract) = {

  show heading: set block(above: 0em, below: 1em)

  let wc = word-count-of(abstract).words

  set align(center)

  linebreak()
  set text(top-edge: 0.7em, bottom-edge: -0.3em)

  v(2em)

  heading("Abstract", outlined: false)

  set align(left)
  set par(justify: true, leading: 1.3em)

  abstract
}
//  // }}}

// init-acknowledgements() // {{{
#let init-acknowledgements(acknowledgements) = {
  show heading: set block(above: 0em, below: 1em)
  heading("Acknowledgements")
  acknowledgements
}
//  // }}}

// init-table-of-contents() // {{{
#let init-table-of-contents() = {
  show heading: set block(above: 0em, below: 0.5em)
  heading("Table of Contents") 
  outline(
    title: []
  )
}
//  // }}}

// init-list-of-tables() // {{{
#let init-list-of-tables() = {
  show heading: set block(above: 0em, below: 0.5em)
  heading("List of Tables")
  outline(
    title: [],
    target: figure.where(kind: table)
  )
}
//  // }}}

// init-list-of-plates() // {{{
#let init-list-of-plates() = {
  show heading: set block(above: 0em, below: 0.5em)
  heading("List of Plates")
  outline(
    title: [],
    target: figure.where(kind: "plate")
  )
}
//  // }}}

// init-list-of-figures() // {{{
#let init-list-of-figures() = {
  show heading: set block(above: 0em, below: 0.5em)
  heading("List of Figures")
  outline(
    title: [],
    target: figure.where(kind: image)
  )
} 
//  // }}}

// init-list-of-appendices() // {{{
#let init-list-of-appendices() = {
  show heading: set block(above: 0em, below: 0.5em)
  heading("List of Appendices")
  outline(
    title: [Appendix],
    target: heading.where(supplement: [Appendix])
  )
}
//  // }}}



// oxthesis() // {{{
#let oxthesis(title: none,
          author: [*missing-param-author*],
          abstract: [],
          acknowledgements: [],
          show-acknowledgements: true,
          show-list-of-tables: true,
          show-list-of-plates: false,
          show-list-of-figures: true,
          show-list-of-appendices: true,
          title-page-top-margin: 5cm,
          title-page-gap-1-height: 4cm,
          title-page-gap-2-height: 4cm,
          title-page-gap-3-height: 5cm,
          title-page-gap-4-height: 3cm,
          title-page-bottom-margin: 3cm,
          page-size-style: "metric",
          main-margin-style: "left-metric",
          font-size: 12pt,
          doc) = {

  let font-size = parse-font(font-size)
  set text(size: font-size)

  show heading: it => block(width: 100%)[
    #set text(weight: "regular")
    #(it)
  ]

  init-title-page(
    title,
    author
  )

  pagebreak()

  set page(numbering: "i")

  init-abstract(
    abstract
  )

  set par(leading: 0.75em)

  pagebreak()

  if (show-acknowledgements) {
    init-acknowledgements(acknowledgements)
    pagebreak()
  }

  init-table-of-contents()

  pagebreak()

  if (show-list-of-tables) {
    init-list-of-tables()
    pagebreak()
  }

  if (show-list-of-plates) {
    init-list-of-plates()
    pagebreak()
  }

  if (show-list-of-figures) {
    init-list-of-figures()
    pagebreak()
  }

  if (show-list-of-appendices) {
    init-list-of-appendices()
    pagebreak()
  }

  set page(numbering: "1")
  counter(page).update(1)

  set heading(numbering: "1.1.1.a")

  doc
}
// // }}}
