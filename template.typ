#import "@preview/physica:0.9.8": innerproduct, evaluated

#let template(it, font-mode: "libertinus") = {
  let body-font = none
  let math-font = none
  let raw-font = none
  let heading-font = none

  if font-mode == "libertinus" {
    body-font = "Libertinus Serif"
    math-font = "Libertinus Math"
    raw-font = "DejaVu Sans Mono"
    heading-font = "Libertinus Sans"
  } else if font-mode == "palatino" {
    body-font = "TeX Gyre Pagella"
    math-font = "TeX Gyre Pagella Math"
    raw-font = "DejaVu Sans Mono"
    heading-font = "TeX Gyre Pagella"
  }
  else if font-mode == "computer" {
    body-font = "New Computer Modern"
    math-font = "New Computer Modern Math"
    raw-font = "New Computer Modern Mono"
    heading-font = "New Computer Modern Sans"
  } else {
    panic("Font mode \"" + font-mode + "\" not valid")
  }

  if font-mode == "computer" {
    show raw: set text(font: raw-font, size: 11pt)
  } else {
    show raw: set text(font: raw-font, size: 9pt)
  }

  set page("us-letter", margin: 1in)

  show math.equation: set text(font: math-font)
  show math.equation: set block(breakable: true)

  show image: set align(center)
  
  // show link: set text(blue)
  // show link: underline

  set par(spacing: 1.2em)

  set text(font: body-font, size: 11pt)
  
  set list(marker: (sym.bullet, sym.bullet.stroked, sym.dash.en))

  set page(numbering: "1")

  show heading: set text(font: heading-font)
  show heading.where(level: 1): it => {
    set text(size: 1.2em)
    block(sticky: true, breakable: false, {
      counter(heading).display(it.numbering)
      h(0.4em)
      it.body 
      v(-0.8em)
      line(length: 100%, stroke: 0.5pt + gray)
      v(0.5em)
    })
  }
  show heading.where(level: 2): set text(size: 1.1em)
  show heading.where(level: 3): set text(size: 1.05em)
  show heading.where(level: 4): it => {
    if font-mode == "libertinus" {
      // Libertinus Sans doesn't have bold italic
      set text(style: "normal", weight: "bold")
    } else {
      set text(style: "italic", weight: "bold")
    }
    
    [
      #counter(heading).display(it.numbering)
      #it.body.
    ]
  }

  let normal-text(it) = {
    set text(font: body-font)
    it
  }
  let nt = normal-text
  
  it
}

#show: template

// Document elements -----------------------
#let title(title, name: "Matthew Ward", date: none, add-divider: true) = {
  set align(center)
  show divider: set line(stroke: 0.5pt + black)

  text(weight: "bold", title)
  
  linebreak()
  text(weight: "medium", name)
  
  linebreak()
  {
    set text(weight: "medium", style: "italic")
    let date-dt = datetime.today()
    if date != none {
      date-dt = toml(bytes("date = " + date)).date
    }
    date-dt.display("[month repr:long] [day padding:none], [year]")
  }
  
  // if add-divider {
  //   block(height: 0.25in, divider())
  // }
}

#let private(it) = {
  set text(blue)
  it
}

#let to-do(it) = {
  set text(fuchsia)
  it
}

#let purple-star = {
  set text(size: 1.5em)
  to-do(sym.star)
}

// recTANGle :)
#let tang(it, color: gray, width: 100%, sticky: false) = {
  align(center, block(
    breakable: true, 
    width: width,
    inset: 1em, 
    stroke: color + 0.5pt,
    sticky: sticky,
    fill: color.lighten(97%),
    {
      set align(left)
      it
    }
  ))
}

#let acknowledgement = tang.with(color: green, width: auto)

#let thread(it, color: luma(50%)) = {
  let border-stroke = stroke(
    paint: color, 
    thickness: 1pt, 
    cap: "square"
  )
  let thread-cap = line(
    stroke: border-stroke, 
    length: 0.25em
  )
  
  block(sticky: true, thread-cap)
  v(0pt, weak: true)
  block(
    breakable: true,
    stroke: (left: border-stroke),
    width: 100%,
    inset: (left: 1em, y: 0.75em),
    spacing: 0pt,
    {
      it
    }
  )
  v(0pt, weak: true)
  block(sticky: false, thread-cap)
}

// Math ----------------------\

#let normal-text(it) = {
  it
}
#let nt = normal-text

#let qed-here = {
  h(1fr)
  sym.square
}

#let proof(it) = [
  _Proof._ #h(0.25em) #it #qed-here
]

#let num-eq(it) = {
  set math.equation(numbering: "(1)")
  it
}

#let exercise(label, it, level: 1) = tang(color: blue, width: 100% - 10% * (level - 1), sticky: true, {
  [==== #label]
  h(0.2em)
  it
})

#let problem = exercise // alias
#let subproblem = exercise.with(level: 2)

#let transpose = $upright(sans(T))$
#let trp = transpose
