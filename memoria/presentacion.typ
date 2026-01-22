#import "@preview/touying:0.6.1": *
#import themes.university: *
#set figure(numbering: none)
#show: university-theme.with(
  aspect-ratio: "16-9",
  config-info(
    title: [Técnicas matemáticas y programación en Python para la resolución de problemas químicos],
    author: [#align(right)[Cristian Angelo Manica Georgiev]],
    date: [#align(left)[#v(0.8em)26-01-2026]],
    institution: [#align(left)[Universitat de Barcelona]],
  ),
)

#title-slide()

== Índice 
#components.adaptive-columns(outline(title: none, indent: 1em))
= Introducción
== Electroquímica
#grid(columns: (1fr, 2fr),[
#image("threeElectrode.jpg", height: 80%)
],[
#image("ButlerVolmerGibbs.png", width: 100%)
])
#grid(columns: (4fr, 3fr, 3fr ),[
  #v(2em)
#image("BVwiki.png", width: 103%)
],[
  #v(2em)
#image("RDEimg.png", width: 100%)
],[
#v(2em)
#image("dischargelot.png", width: 120%, height: 60%)
]
)
#grid(columns: (1fr, 1fr),[
#image("CVwiki.png", width: 100%)
],[
#image("CVplots.png", height: 105%)
])


== Métodos numéricos
#grid(columns: (1fr, 1fr),[
#image("EulerMethod.png", height: 70%)
],[
#v(1em)
#align(right)[
#block(width: auto) 
#align(center)[Tipos:]
#align(left)[#h(3em)$ cases(
"Integración numérica",
"Derivación numérica",
"Álgebra lineal numérica",
"Resolucionadores de EDOs",
"Resolucionadores de EDPs", 
"Interpolación y ajuste de curvas", 
)$]
]]

)
#grid(columns: (1fr, 1fr),[
#set text(18pt)
#figure(
  block(
    [
      \
      \
      $ x_(n+1) = x_n - f(x_n) / (f'(x_n)) $
    ]
  ),
  caption: "Newton-Raphson"
)
],[
#set text(18pt)
#figure(
  block(
    [
      \
$ mat(b_1, c_1, , , 0; a_2, b_2, c_2, ; , a_3, b_3, dots.down; 
, , dots.down, dots.down, c_(n-1); 0, , , a_n, b_n) vec(x_1, x_2, x_3, dots.v, x_n) = 
vec(d_1, d_2, d_3, dots.v, d_n) $
    ]
  ),
  caption: "Algoritmo de Thomas"
)
],[
#set text(18pt)
#figure(
  block(
    [

\ $cases(
  k_1 = f(x_n, y_n),
  k_2 = f(x_n + (Delta x)/2, y_n + (Delta x)/2 k_1),
  k_3 = f(x_n + (Delta x)/2, y_n + (Delta x)/2 k_2),
  k_4 = f(x_n + Delta x, y_n + Delta x k_3)) $
$ y_(n+1) = y_n + (Delta x)/6 (k_1 + 2k_2 + 2k_3 + k_4) $
    ]
  ),
  caption: "Runge-Kutta 4"
)
],[
#set text(18pt)
#figure(
  block(
    [
      \
      \
      \
$ f''(x) approx (f(x + Delta x) - 2f(x) + f(x - Delta x)) / (Delta x^2) $
    ]
  ),
  caption: "Fórmula de tres puntos"
)
])

= Objetivos
#block(width: auto)[
#set text(size: 24pt)
- Enfoque Pedagógico y Didáctico
- Desarrollo técnico de sistemas electroquímicos
- Consolidación del Conocimiento]

= Resultados
== Implementación
#grid(columns: (1fr, 1fr, 1fr), gutter: 20pt,
[
#set text(size: 18pt)
#figure(
image("ButlerVolmerStructure.png", width: 100%),
caption: "ButlerVolmer.py"
  )
],[
#set text(size: 18pt)
#figure(
image("RotatingDiskElectrodeStructure.png", width: 100%),
caption: "RotatingDiskElectrode.py"
  )
],[
#set text(size: 18pt)
  #figure(
image("BatteryDischargeStructure.png", width: 100%),
    caption: "BatteryDischarge.py"
  )
]
)
#figure(
image("CyclicVoltammetryStructure.png", height: 70%),
    caption: "CyclicVoltammetry.py"
    )

== ButlerVolmer.py
#grid(columns: (1fr, 1fr),[
#set text(size: 17pt)
#figure(
image("test1.png", width: 75%),
caption: [$\ "Fc"^+ + e^- -> "Fc" "with" \ ["Fc"^+] =  10^(-4) "mol/"m^3, \ ["Fc"] = 10^(-4) "mol/"m^3$]
)
],[
#set text(size: 17pt)
#figure(
image("test2.png", width: 75%),
caption: [$\ "Fc"^+ + e^- -> "Fc" "with" \ ["Fc"^+] = 10^(-3) "mol/"m^3, \ ["Fc"] = 10^(-4) "mol/"m^3$]
)
])
#v(-0.5em)
#align(center)[
  #set text(size:20pt)
  #block(
    fill: rgb("#e1f5fe"), 
    stroke: 1pt + blue,
    inset: 12pt,
    radius: 4pt,
    width: 60%,
    [
    $ j = j_0 (exp(((1 - beta) n F eta)/(R T)) - exp(-(beta n F eta)/(R T))) $
    ]
  )]
== RotatingDiskElectrode.py
#v(1em)
#grid(columns: (1fr, 1fr, 1fr), [
#set text(size: 18pt)
#figure(
image("rpmtest2.png", width: 95%),
caption: [Standard limit intensity measurement with ferrocene \ $i_L = 2 · 10^(-4)A$]
)
],[
#set text(size: 18pt)
#figure(
image("rpmtest3.png", width: 95%),
caption: [Limit current for 10 times the concentration of reactants \ $i_L = 2 · 10^(-3)A$]
)
],[
#set text(size: 18pt)
#figure(
image("rpmtest4.png", width: 95%),
caption: [Limit current for 10 times the area of the electrode \ $i_L = 2 · 10^(-3)A$]
)
]
)
#align(center)[#block(
    fill: rgb("#e1f5fe"), 
    stroke: 1pt + blue,
    inset: 12pt,
    radius: 4pt,
    width: 45%,
    [
      $ i_L = 0.620n F A D^(2/3) omega^(1/2) nu^(-1/6) C $
    ]
  )]
== BatteryDischarge.py
#grid(columns: (1fr, 1fr, 1fr),[
#set text(size: 18pt)
#figure(
image("battest1.png", width: 95%),
caption: [Battery Discharge]
)
],[
#set text(size: 18pt)
#figure(
image("battest2.png", width: 95%),
caption: [Battery Discharge for 10 $I R_"int"$]
)
],[
#set text(size: 18pt)
#figure(
image("battest3.png", width: 95%),
caption: [Battery Discharge for $"1/5" Q_"nom"$]
)
])
#grid(columns: (1fr, 1fr),[
  #align(center)[#block(
    fill: rgb("#e1f5fe"), 
    stroke: 1pt + blue,
    inset: 32pt,
    radius: 4pt,
    width: 75%,
    [
$ V = "OCV"("SOC") - I R_"int" $
    ]
  )
  ]],[
  #align(center)[#block(
    fill: rgb("#e1f5fe"), 
    stroke: 1pt + blue,
    inset: 12pt,
    radius: 4pt,
    width: 110%,
    [
    $ "SOC"(t) = "SOC"_(t = 0) - 1 / Q_"nom"integral^t_0 I(t) d t $
    ]
  )
]])
== CyclicVoltammetry.py
#grid(columns: (1fr, 1fr),[
#set text(size: 18pt)
#figure(
image("CVtest1.png", width: 50%),
caption: [Voltammogram, ferrocene, \ $C_(O,"bulk") = 1 "mM", C_(R,"bulk") = 1 "mM"$]
)
],[
#set text(size: 18pt)
#figure(
image("CVtest2.png", width: 50%),
caption: [Voltammogram, $1000 k_0$ ]
)
],[
  #v(2em)
  #align(center)[#block(
    fill: rgb("#e1f5fe"), 
    stroke: 1pt + blue,
    inset: 12pt,
    radius: 4pt,
    width: 50%,
    [
$ J = - D (partial C) / (partial x) $
    ])]
],[
  #v(2em)
  #align(center)[#block(
    fill: rgb("#e1f5fe"), 
    stroke: 1pt + blue,
    inset: 12pt,
    radius: 4pt,
    width: 50%,
    [
$ (partial C) / (partial t) = D (partial^2 C) / (partial x^2) $
    ])]
  ],[
#set text(size: 18pt)
#figure(
image("CVtest3.png", width: 50%),
caption: [Voltammogram, $beta_c = 0.75$]
)
],[
#set text(size: 18pt)
#figure(
image("CVtest4.png", width: 50%),
caption: [Voltammogram, $E_"start"$ not at the equilibrium potential]
)
],[
  #v(2em)
  #align(center)[#block(
    fill: rgb("#e1f5fe"), 
    stroke: 1pt + blue,
    inset: 12pt,
    radius: 4pt,
    width: 50%,
    [
$ J = - D (partial C) / (partial x) $
    ])]
],[
  #v(2em)
  #align(center)[#block(
    fill: rgb("#e1f5fe"), 
    stroke: 1pt + blue,
    inset: 12pt,
    radius: 4pt,
    width: 50%,
    [
$ (partial C) / (partial t) = D (partial^2 C) / (partial x^2) $
    ])]
  ],[
#set text(size: 18pt)
#figure(
image("CVtest3.png", width: 50%),
caption: [Voltammogram, $beta_c = 0.75$]
)
],[
#set text(size: 18pt)
#figure(
image("CVtest4.png", width: 50%),
caption: [Voltammogram, $E_"start"$ not at the equilibrium potential]
)
],[
  #v(2em)
  #align(center)[#block(
    fill: rgb("#e1f5fe"), 
    stroke: 1pt + blue,
    inset: 12pt,
    radius: 4pt,
    width: 50%,
    [
$ J = - D (partial C) / (partial x) $
    ])]
],[
  #v(2em)
  #align(center)[#block(
    fill: rgb("#e1f5fe"), 
    stroke: 1pt + blue,
    inset: 12pt,
    radius: 4pt,
    width: 50%,
    [
$ (partial C) / (partial t) = D (partial^2 C) / (partial x^2) $
    ])]
  ],[
#set text(size: 18pt)
#figure(
image("CVtest3.png", width: 50%),
caption: [Voltammogram, $beta_c = 0.75$]
)
],[
#set text(size: 18pt)
#figure(
image("CVtest4.png", width: 50%),
caption: [Voltammogram, $E_"start"$ not at the equilibrium potential]
)
],[
  #v(2em)
  #align(center)[#block(
    fill: rgb("#e1f5fe"), 
    stroke: 1pt + blue,
    inset: 12pt,
    radius: 4pt,
    width: 50%,
    [
$ J = - D (partial C) / (partial x) $
    ])]
],[
  #v(2em)
  #align(center)[#block(
    fill: rgb("#e1f5fe"), 
    stroke: 1pt + blue,
    inset: 12pt,
    radius: 4pt,
    width: 50%,
    [
$ (partial C) / (partial t) = D (partial^2 C) / (partial x^2) $
    ])]],[
#set text(size: 18pt)
#figure(
image("pres1.png", width: 100%),
caption: [$t = t_0$]
)
],[
#set text(size: 18pt)
#figure(
image("pres2.png", width: 100%),
caption: [Voltammogram, $E_"start"$ not at the equilibrium potential]
)
  ],[
#set text(size: 18pt)
#figure(
image("pres3.png", width: 100%),
caption: [Voltammogram, $beta_c = 0.75$]
)
],[
#set text(size: 18pt)
#figure(
image("pres4.png", width: 100%),
caption: [Voltammogram, $E_"start"$ not at the equilibrium potential]
)
  ],[
#set text(size: 18pt)
#figure(
image("pres5.png", width: 100%),
caption: [Voltammogram, $beta_c = 0.75$]
)
],[
#set text(size: 18pt)
#figure(
image("pres6.png", width: 100%),
caption: [Voltammogram, $E_"start"$ not at the equilibrium potential]
)
],[
#set text(size: 18pt)
#figure(
image("pres7.png", width: 100%),
caption: [Voltammogram, $E_"start"$ not at the equilibrium potential]
)
],[
#set text(size: 18pt)
#figure(
image("pres8.png", width: 100%),
caption: [Voltammogram, $E_"start"$ not at the equilibrium potential]
)
  ])

= Conclusions
#block(width: auto)[
#set text(size: 24pt)
- Development of functional software
- Precision through numerical methods
- Achievement of the pedagogical objective]
#pagebreak()
#align(center + horizon)[#set text(size: 70pt); *GRACIAS POR  SU ATENCIÓN*]
#grid(columns: (1fr, 1fr), gutter: 8pt, [
#image("BV1.png", width: 100%)
],[
#image("BV2.png", width: 100%)
],[
#image("RotatingDiskElectrode.png", width: 100%)
],[
#image("BatteryDischarge.png", height: 100%)
],[
#image("CVp2.png", width: 100%)
],[
#image("CVp1.png", width: 100%)
],[
#image("CVp3.png", width: 100%)
],[
#image("CVp4.png", width: 100%)
],[
#image("CVp5.png", width: 100%)
],[
#image("CVp6.png", width: 100%)
],[
#image("CVp7.png", height: 100%)
],[
#image("CVp8.png", width: 100%)
])
