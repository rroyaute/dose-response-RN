// Some definitions presupposed by pandoc's typst output.
#let blockquote(body) = [
  #set text( size: 0.92em )
  #block(inset: (left: 1.5em, top: 0.2em, bottom: 0.2em))[#body]
]

#let horizontalrule = [
  #line(start: (25%,0%), end: (75%,0%))
]

#let endnote(num, contents) = [
  #stack(dir: ltr, spacing: 3pt, super[#num], contents)
]

#show terms: it => {
  it.children
    .map(child => [
      #strong[#child.term]
      #block(inset: (left: 1.5em, top: -0.4em))[#child.description]
      ])
    .join()
}

// Some quarto-specific definitions.

#show raw.where(block: true): block.with(
    fill: luma(230), 
    width: 100%, 
    inset: 8pt, 
    radius: 2pt
  )

#let block_with_new_content(old_block, new_content) = {
  let d = (:)
  let fields = old_block.fields()
  fields.remove("body")
  if fields.at("below", default: none) != none {
    // TODO: this is a hack because below is a "synthesized element"
    // according to the experts in the typst discord...
    fields.below = fields.below.amount
  }
  return block.with(..fields)(new_content)
}

#let empty(v) = {
  if type(v) == "string" {
    // two dollar signs here because we're technically inside
    // a Pandoc template :grimace:
    v.matches(regex("^\\s*$")).at(0, default: none) != none
  } else if type(v) == "content" {
    if v.at("text", default: none) != none {
      return empty(v.text)
    }
    for child in v.at("children", default: ()) {
      if not empty(child) {
        return false
      }
    }
    return true
  }

}

// Subfloats
// This is a technique that we adapted from https://github.com/tingerrr/subpar/
#let quartosubfloatcounter = counter("quartosubfloatcounter")

#let quarto_super(
  kind: str,
  caption: none,
  label: none,
  supplement: str,
  position: none,
  subrefnumbering: "1a",
  subcapnumbering: "(a)",
  body,
) = {
  context {
    let figcounter = counter(figure.where(kind: kind))
    let n-super = figcounter.get().first() + 1
    set figure.caption(position: position)
    [#figure(
      kind: kind,
      supplement: supplement,
      caption: caption,
      {
        show figure.where(kind: kind): set figure(numbering: _ => numbering(subrefnumbering, n-super, quartosubfloatcounter.get().first() + 1))
        show figure.where(kind: kind): set figure.caption(position: position)

        show figure: it => {
          let num = numbering(subcapnumbering, n-super, quartosubfloatcounter.get().first() + 1)
          show figure.caption: it => {
            num.slice(2) // I don't understand why the numbering contains output that it really shouldn't, but this fixes it shrug?
            [ ]
            it.body
          }

          quartosubfloatcounter.step()
          it
          counter(figure.where(kind: it.kind)).update(n => n - 1)
        }

        quartosubfloatcounter.update(0)
        body
      }
    )#label]
  }
}

// callout rendering
// this is a figure show rule because callouts are crossreferenceable
#show figure: it => {
  if type(it.kind) != "string" {
    return it
  }
  let kind_match = it.kind.matches(regex("^quarto-callout-(.*)")).at(0, default: none)
  if kind_match == none {
    return it
  }
  let kind = kind_match.captures.at(0, default: "other")
  kind = upper(kind.first()) + kind.slice(1)
  // now we pull apart the callout and reassemble it with the crossref name and counter

  // when we cleanup pandoc's emitted code to avoid spaces this will have to change
  let old_callout = it.body.children.at(1).body.children.at(1)
  let old_title_block = old_callout.body.children.at(0)
  let old_title = old_title_block.body.body.children.at(2)

  // TODO use custom separator if available
  let new_title = if empty(old_title) {
    [#kind #it.counter.display()]
  } else {
    [#kind #it.counter.display(): #old_title]
  }

  let new_title_block = block_with_new_content(
    old_title_block, 
    block_with_new_content(
      old_title_block.body, 
      old_title_block.body.body.children.at(0) +
      old_title_block.body.body.children.at(1) +
      new_title))

  block_with_new_content(old_callout,
    new_title_block +
    old_callout.body.children.at(1))
}

// 2023-10-09: #fa-icon("fa-info") is not working, so we'll eval "#fa-info()" instead
#let callout(body: [], title: "Callout", background_color: rgb("#dddddd"), icon: none, icon_color: black) = {
  block(
    breakable: false, 
    fill: background_color, 
    stroke: (paint: icon_color, thickness: 0.5pt, cap: "round"), 
    width: 100%, 
    radius: 2pt,
    block(
      inset: 1pt,
      width: 100%, 
      below: 0pt, 
      block(
        fill: background_color, 
        width: 100%, 
        inset: 8pt)[#text(icon_color, weight: 900)[#icon] #title]) +
      if(body != []){
        block(
          inset: 1pt, 
          width: 100%, 
          block(fill: white, width: 100%, inset: 8pt, body))
      }
    )
}



#let article(
  title: none,
  authors: none,
  date: none,
  abstract: none,
  abstract-title: none,
  cols: 1,
  margin: (x: 1.25in, y: 1.25in),
  paper: "us-letter",
  lang: "en",
  region: "US",
  font: (),
  fontsize: 11pt,
  sectionnumbering: none,
  toc: false,
  toc_title: none,
  toc_depth: none,
  toc_indent: 1.5em,
  doc,
) = {
  set page(
    paper: paper,
    margin: margin,
    numbering: "1",
  )
  set par(justify: true)
  set text(lang: lang,
           region: region,
           font: font,
           size: fontsize)
  set heading(numbering: sectionnumbering)

  if title != none {
    align(center)[#block(inset: 2em)[
      #text(weight: "bold", size: 1.5em)[#title]
    ]]
  }

  if authors != none {
    let count = authors.len()
    let ncols = calc.min(count, 3)
    grid(
      columns: (1fr,) * ncols,
      row-gutter: 1.5em,
      ..authors.map(author =>
          align(center)[
            #author.name \
            #author.affiliation \
            #author.email
          ]
      )
    )
  }

  if date != none {
    align(center)[#block(inset: 1em)[
      #date
    ]]
  }

  if abstract != none {
    block(inset: 2em)[
    #text(weight: "semibold")[#abstract-title] #h(1em) #abstract
    ]
  }

  if toc {
    let title = if toc_title == none {
      auto
    } else {
      toc_title
    }
    block(above: 0em, below: 2em)[
    #outline(
      title: toc_title,
      depth: toc_depth,
      indent: toc_indent
    );
    ]
  }

  if cols == 1 {
    doc
  } else {
    columns(cols, doc)
  }
}

#set table(
  inset: 6pt,
  stroke: none
)
#show: doc => article(
  title: [Dose-response models as individual reaction norms: Acocunting for individual variation in ecotoxicological assays],
  date: [Tuesday, August 27, 2024],
  toc: true,
  toc_title: [Contents],
  toc_depth: 3,
  cols: 1,
  doc,
)


== Problem description and model specification
<problem-description-and-model-specification>
Phenotypic variation is common in biology, yet is poorly accounted for in ecotoxicology and environmental risk assessment. When variation is ignored, this may lead to decision that do not protect the most sensitive individuals in a population. The following study attempts to answer the following fundamental questions:

+ What is the risk of ignoring biological variation when estimating ecotoxicological parameters?
+ What are optimal designs and sample sizes for estimating among-individual or among-genotype variation?

To answer (1), we run a simulation study where we compare the estimation of ecotoxicological parameters while varying the magnitude of among-genotype differences. To answer (2), we simulate data with increased sample sizes and / or different experimental designs to determine which design yields the best estimation of among-individual differences in ecotoxicological parameters

== A hierarchical dose-response model for capturing heterogeneous responses to contaminant exposure
<a-hierarchical-dose-response-model-for-capturing-heterogeneous-responses-to-contaminant-exposure>
We chose a common dose-response model where a biological response $R$ remains unchanged while the contaminant dose remains below a No Effect Concentration (NEC). As soon as the dose falls above the NEC, the response declines at an exponential rate $- e^beta$ until it reaches a lower bound $R_(m i n)$.

$ R = R_(m i n) + (R_(m a x) - R_(m i n)) e^(- e^(beta times (D o s e - N E C) times (D o s e > N E C))) $ #box(image("../outputs/figs/fig_DR_NEC.png"))

We assume the response is positive with an average value $R_(m a x) = 100$ and a lower bound $R_(m i n) = 1$. We then assume that we are able to observe the response of distinct genotypes, each varying in their sensitivity to the contaminant. This means that each genotype can have a unique value for the parameters $R_(m a x)$, $N E C$ and $beta$. For simplicity, we assume that the lower bound $R_(m i n)$ is shared among all genotypes.

The model can therefore be rewritten as a hierarchical model such that the average value of each genotype differs from the population average by a normally-distributed offset that depends on the amount of among-genotype variation. We use a lognormal likelihood to contrain the estimations to be stricly positive and use a Bayesian framework with weak priors to constrain the data around plausible values. The full model can thus be written as:

\$\$
\$\$

Here is what the patterns would look like assuming an among-genotype sd of 0.5:

#figure([
#box(image("../outputs/figs/fig.param.sim.3.jpeg"))
], caption: figure.caption(
position: bottom, 
[
Fig 2. Dose-response with 1000 genotypes deviating from the average trend for a given parameters
]), 
kind: "quarto-float-fig", 
supplement: "Figure", 
)


== When is individual variation problematic for the estimation of ecotoxicological parameters?
<when-is-individual-variation-problematic-for-the-estimation-of-ecotoxicological-parameters>
With this general model, we will investigate three cases for which the presence of among-genotype variance can potentially bias our estimation of the population averages for the ecotoxicological parameters $mu_(N E C)$ and $mu_beta$:

- Case 1: Only the average response $R_(m a x)$ differs among genotypes
- Case 2: Only the $N E C$ differs among genotypes
- Case 3: Only the decay rate $beta$ differs among genotypes
- Case 4: All three parameters vary independently among genotypes
- Case 5: All three parameters covary among genotypes

In each case, we generate 1000 datasets by simulating from the model version including the among-genotype variation for the specific set of parameters under consideration. We compare the fit of a model ignoring individual difference with a correctly-specified model. We rely on the following metrics to assess our models performances:

- % Bias: a metric indicating how far off our estimate is from the true parameter value and defined as the percent distance to the true value of the parameter $(theta - hat(theta)) / hat(theta) times 100$
- Precision: a metric indicating the agreement in parameter values for a given set of simulations and defined as the difference between the 25th and 75th quantiles for a given set of simulations $q_(s i m_25) - q_(s i m_75)$
- RMSE: an additional metric for the distance between the true value and the model estimate, defined as the square root of average difference between the true value of the parameter and the posterior median estimated from the model $sqrt(E (theta - hat(theta)))$

We also want to understand how strong these genotype differences need to be so that the estimation of ecotoxicological parameters is problematic when not explicitly modeled. We will therefore run each case with increasing values for the among-genotype parameters: $C V_i = [0 , 1 , 10 , 20 , 30 , 40 , 50]$, where $C V_i$ corresponds to the among-genotype coefficient of variation, calculated as $sigma_(i_theta) = C V_(i_theta) times mu_theta$.

The full simulation plan is summarized in the table below:

#show figure: set block(breakable: true)

#let nhead = 1;
#let nrow = 4;
#let ncol = 4;

  #let fill-array = ( 
    // tinytable cell fill after
    (y: 0, x: 3, fill: black),
    (y: 0, x: 2, fill: black),
    (y: 0, x: 1, fill: black),
    (y: 0, x: 0, fill: black),
  )
  #let style-array = ( 
    // tinytable cell style after
    (y: 0, x: 3, color: white, underline: false, italic: false, bold: false, mono: false, strikeout: false, fontsize: 1em, indent: false),
    (y: 0, x: 2, color: white, underline: false, italic: false, bold: false, mono: false, strikeout: false, fontsize: 1em, indent: false),
    (y: 0, x: 1, color: white, underline: false, italic: false, bold: false, mono: false, strikeout: false, fontsize: 1em, indent: false),
    (y: 0, x: 0, color: white, underline: false, italic: false, bold: false, mono: false, strikeout: false, fontsize: 1em, indent: false),
  )
  #let align-array = (
    // tinytable cell align after
    (y: 0, x: 3, align: center),
    (y: 0, x: 2, align: center),
    (y: 0, x: 1, align: center),
    (y: 0, x: 0, align: center),
    (y: 4, x: 3, align: center),
    (y: 3, x: 3, align: center),
    (y: 2, x: 3, align: center),
    (y: 1, x: 3, align: center),
    (y: 0, x: 3, align: center),
    (y: 4, x: 2, align: center),
    (y: 3, x: 2, align: center),
    (y: 2, x: 2, align: center),
    (y: 1, x: 2, align: center),
    (y: 0, x: 2, align: center),
    (y: 4, x: 1, align: center),
    (y: 3, x: 1, align: center),
    (y: 2, x: 1, align: center),
    (y: 1, x: 1, align: center),
    (y: 0, x: 1, align: center),
    (y: 4, x: 0, align: center),
    (y: 3, x: 0, align: center),
    (y: 2, x: 0, align: center),
    (y: 1, x: 0, align: center),
    (y: 0, x: 0, align: center),
  )
  // tinytable align-default-array after
  #let align-default-array = ( left, left, left, left, ) // tinytable align-default-array here
  #show table.cell: it => {
    let tmp = it
    let data = style-array.find(data => data.x == it.x and data.y == it.y)
    if data != none {
      set text(data.color)
      set text(data.fontsize)
      if data.indent != false { tmp = pad(left: data.indent, tmp) }
      if data.underline == true { tmp = underline(tmp) }
      if data.italic == true { tmp = emph(tmp) }
      if data.bold == true { tmp = strong(tmp) }
      if data.mono == true { tmp = math.mono(tmp) }
      if data.strikeout == true { tmp = strike(tmp) }
      tmp
    } else {
      tmp
    }
  }

  #align(center, [

  #table( // tinytable table start
    columns: (auto, auto, auto, auto),
    stroke: none,
    align: (x, y) => {
      let data = align-array.find(data => data.x == x and data.y == y)
      if data != none {
        data.align
      } else {
        align-default-array.at(x)
      }
    },
    fill: (x, y) => {
      let data = fill-array.find(data => data.x == x and data.y == y)
      if data != none {
        data.fill
      }
    },

    // tinytable lines after
table.hline(y: 5, start: 0, end: 4, stroke: 0.1em + black),
table.hline(y: 1, start: 0, end: 4, stroke: 0.05em + black),
table.hline(y: 0, start: 0, end: 4, stroke: 0.1em + black),

    table.header(
      repeat: true,
[Parameter], [mean], [sd], [Simulations],
    ),

    // tinytable cell content after
[$R_{min}$], [  1], [$0$                                                    ], [1000],
[$R_{max}$], [100], [$[0, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5]$                   ], [1000],
[$NEC$    ], [ 10], [$[0, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5] \times \mu_{NEC}$], [1000],
[$\beta$ ], [  0.01], [$[0, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5]$                   ], [1000],

    table.footer(
      repeat: false,
      // tinytable notes after
    ),

  ) // end table

  ]) // end align
== Limits of the dose-response reaction norm approach
<limits-of-the-dose-response-reaction-norm-approach>
In many cases, we do not have access to genotypes or clonal lines we can neatly expose through the totality of the dose gradient. What is much more common is to expose sample groups of individuals and randomly assign them to a dose treatment. In this section we consider three additional scenarios that more closely match with the way ecotoxicological data is gathered:

+ Individual responses are measured only once and for a unique dose along the gradient
+ Individual responses are measured multiple times during the course of exposure to a single dose
+ Individual responses are measured repeatedly before and after exposure

For each scenario we ask the following questions: (#emph[i];) What are optimal sampling design for a robust estimation of ecotoxicological parameters ($N E C$ and $beta$) and (#emph[ii];) How should individual differences in be accounted for?

=== What are optimal sampling designs to deal with individual differences in ecotoxicology in absence of repeated measurements?
<what-are-optimal-sampling-designs-to-deal-with-individual-differences-in-ecotoxicology-in-absence-of-repeated-measurements>
=== What are optimal sampling designs to estimate individual differences in ecotoxicology?
<what-are-optimal-sampling-designs-to-estimate-individual-differences-in-ecotoxicology>



