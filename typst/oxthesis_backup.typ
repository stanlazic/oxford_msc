#import "@preview/wordometer:0.1.4": word-count, total-words, word-count-of, string-word-count
#show: word-count

// parse-page-dims() // {{{
#let parse-page-dims(page-size-style) = {
  let (page-width, page-height) = if page-size-style == "metric" {
    (21.5cm, 28cm)
  } else if (page-size-style == "imperial") {
    (8.5in, 11in)
  } else {
    panic([Parameter `page-size-style` must be either `metric` or `imperial`])
  }
  (page-width, page-height)
}
//  // }}}

// parse-margin-dims() // {{{
#let parse-margin-dims(page-margin-style) = {
  let (margin-t, margin-b, margin-l, margin-r) = if (page-margin-style == "left-metric") {
    (20mm, 20mm, 32mm, 20mm)
  } else if (page-margin-style == "left-imperial") {
    (0.75in, 0.75in, 1.25in, 0.75in)
  } else if (page-margin-style == "metric") {
    (20mm, 20mm, 20mm, 20mm)
  } else if (page-margin-style == "imperial") {
    (0.75in, 0.75in, 0.75in, 0.75in)
  } else {
    panic([Parameter `page-margin-style` must be one of `left-metric`, `left-imperial`, `metric`, or `imperial`])
  }
  (margin-t, margin-b, margin-l, margin-r)
}

//  // }}}

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
                author,
                title-page-top-margin,
                title-page-gap-1-height,
                title-page-gap-2-height,
                title-page-gap-3-height,
                title-page-gap-4-height,
                title-page-bottom-margin,
                page-size-style
              ) = {
  set par(spacing: 0em)

  let top-margin = if title-page-top-margin in (2in, 5cm) {
    title-page-top-margin
  } else {
    panic([Parameter `title-page-top-margin` must be either `2in` or `5cm`])
  }

  let gap-1-height = if title-page-gap-1-height in (1.5in, 4cm) {
    title-page-gap-1-height
  } else {
    panic([Parameter `title-page-gap-1-height` must be either `1.5in` or `4cm`])
  }

  let gap-2-height = if title-page-gap-2-height in (1.5in, 4cm) {
    title-page-gap-2-height
  } else {
    panic([Parameter `title-page-gap-2-height` must be either `1.5in` or `4cm`])
  }

  let gap-3-height = if title-page-gap-3-height in (2in, 5cm) {
    title-page-gap-3-height
  } else {
    panic([Parameter `title-page-gap-3-height` must be either `2in` or `5cm`])
  }

  let gap-4-height = if title-page-gap-4-height in (1.25in, 3cm) {
    title-page-gap-4-height
  } else {
    panic([Parameter `title-page-gap-4-height` must be either `1.25in` or `3cm`])
  }

  let bottom-margin = if title-page-bottom-margin in (1.25in, 3cm) {
    title-page-bottom-margin
  } else {
    panic([Parameter `title-page-bottom-margin` must be either `1.25in` or `3cm`])
  }

  let (width, height) = parse-page-dims(page-size-style)
  set page(
    page: "a4",		
    width: width,    
    height: height,
    margin: (
      top: top-margin,
      bottom: bottom-margin,
      left: 32mm,
      right: 32mm
    ),
  ) 

  align(center)[
    #text(size: 1.4em)[#title]
    #v(gap-1-height)
    #image("img/beltcrest.svg")
    #v(gap-2-height)
    #author
    #v(gap-1-height)
    This dissertation has been submitted to the University of Oxford in partial fulfilment of the
requirement for the award of the degree of MSc in Evidence-Based Health Care (Medical
Statistics).
    #v(gap-1-height)
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
  let (page-width, page-height) = parse-page-dims(
    page-size-style
  )
  let (margin-t, margin-b, margin-l, margin-r) = parse-margin-dims(
    main-margin-style
  )
  set page(
    width: page-width,
    height: page-height,
    margin: (top: margin-t, bottom: margin-b, left: margin-l, right: margin-r),
  )
  let font-size = parse-font(font-size)
  set text(size: font-size)

  show heading: it => block(width: 100%)[
    #set text(weight: "regular")
    #(it)
  ]

  init-title-page(
    title,
    author,
    title-page-top-margin,
    title-page-gap-1-height,
    title-page-gap-2-height,
    title-page-gap-3-height,
    title-page-gap-4-height,
    title-page-bottom-margin,
    page-size-style,
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
