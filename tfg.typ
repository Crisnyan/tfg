#set document(title: [Manica_Georigiev_TFG])
#let title = "Mathematical techniques and Python programming for solving chemical problems."
#set page(
  paper: "a4",
  margin: (top: 2cm, bottom: 2cm, left: 60pt, right: 2.5cm),
)

#set text(
  font: "Arial",
  size: 10pt,
  stretch: 80%
)
#set par(
  justify: true
)

#let weird(big, small, t_big, t_small) = {
  text(size: big)[*#t_big*]
  text(size: small)[*#t_small*]
}

// NOTE: portada
#box(width: auto)[
  #align(right)[#text(fill: gray.darken(50%), size: 9pt)[Tutor/s]] 
  #v(-0.5em)
#line(length: 30%, stroke: 0.5pt + gray.darken(50%))
  #v(-0.5em)
  #align(right)[Dr. Miguel Gonzalez Perez]
  #v(-0.5em)
  #align(right)[#text(style: "italic", size: 9pt)[Departament de Química Xxxxx]]
  #v(0.2em)
]

#v(0.5cm)
#align(left)[
  #image("grauUB.png", width: 5cm)
]

#v(1.6cm)


#align(left)[
  #block(width: 100%)[
    #text(weight: "bold", size: 10pt)[
    #title \
    Técnicas matemáticas y programación en Python para la resolución de problemas químicos.
]]]

#v(1.5cm)

#align(center)[
  #text(size: 24pt, weight: 800, stretch: 100%)[*Treball Final de Grau*]
]

#align(left)[
  #text(size: 14pt)[Cristian Manica Georgiev] \
  #h(2cm)
  #text(size: 10pt)[January 2026]
]

#align(center)[#image("logoUB.png", width: 40%)]

#v(2cm)

#align(center)[#box(width: auto)[
#align(left)[
    Aquesta obra està subjecta a la llicència de: \
    Reconeixement NoComercial-SenseObraDerivada
  #v(-0.75em)
  #image("cc.png", width: 2.75cm)
  #v(-0.75em)
  #text(size: 10pt)[#link("http://creativecommons.org/licenses/by-nc-nd/3.0/es/")]
]]]

#pagebreak()

// NOTE: agradecimientos
Quiero agradecer primero a mi madre y a mi padre, por hacer posible mis estudios universitarios, 
respetando mi decisión de estudiar química. \
Agradezco todo el apoyo de todas las personas que han hecho posible que llegue a ser la persona 
que soy hoy: mis amigos, ayudándome cuando tengo problemas, los profesores de la facultad por 
enseñarme la química y demostrarme que no se límita a la triste química que me enseñaron en 
bachillerato, y por último, a mi tutor por ayudarme con este TFG, muchas gracias por toda la 
ayuda otorgada a todos, no seria nadie sin vosotros.

#pagebreak()
// NOTE: Report

#align(center)[#h(4cm)#set text(fill: black.lighten(35%))
#weird(36pt, 26pt, "R", "EPORT")]

#pagebreak()

// NOTE: SDG

#weird(16pt, 11pt, "I", "DENTIFICATION AND REFLECTION ON THE SUSTAINABLE DEVELOPMENT GOALS") 
#weird(16pt, 11pt, "(SDG)", "") \

#v(1em)
The main objective of this TFG relates to the 5 Ps of the Sustainable Development 
Goals in two ways:
Fist, it relates to people by providing computational tools, thus facilitating
high-level scientific and technical training in electrochemistry, aiding students
with the visualization of the theoretical equations, providing a visual representation
that furthers their understainding on the math.
And second, it focuses on the planet by providing alternatives to chemical experiments, 
making possible the prevention of generating chemical waste, which could be hazarduous. \
This work aligns with three specific SDG goals and targets: SDG 4, SDG 9 and SDG 12. \
Quality education (SDG 4) is reflected in the before mentioned P of the SDGs, People.
The software supports target 4.4 by increasing the number of people who have vocation
for chemistry, allowing them to further their skills for employment in the chemical and 
digital industries. \
Industry, innovation and infrastructure (SDG 5) is shown in this work by providing 
a simple numerical solver, using RK4 and implicit methods, which as point 9.5 implies, 
it's a way of enhancing industrical technological capabilities. \
To finalize, respondible consumption and production (SDG 12) is apparent in this work 
by its adherence to the point 12.4, adhering to green chemistry principles which reduce 
the release of chemicals into the air, water and other mediums. \
Should this project grow beyond its scope, the final projection would be the optimization 
and adoption into an open-access educational platform, which would benefit enormously 
students who have a hard time grasping the meaning behind the mathematical equations in 
chemistry. 

#pagebreak()

// NOTE: Contents

#let page_number = 1

#let which_header(n) = {
  if calc.odd(n) {
    [#emph[#text(size: 7pt, stretch: 70%)[#title] #n] #v(-0.5em)]
    v(-0.4em)
    line(length: 100%, stroke: 0.5pt + gray.darken(50%))
  } else {
    [#emph[#n #text(size: 7pt, stretch: 70%)[Manica Georgiev, Cristian Angelo]] #v(-0.3em)]
    v(-0.4em)
    line(length: 100%, stroke: 0.5pt + gray.darken(50%))
  }
}

#set page(header: which_header(page_number))
#set text(stretch: 100%)
#v(2.5cm)
#weird(16pt, 11pt, "C", "ONTENTS")





#pagebreak()

// NOTE: Página en blanco

#let page_number = page_number + 1
#set page(header: which_header(page_number))
#page([])
#pagebreak()

// NOTE: Summary

#let page_number = page_number + 1
#set page(header: which_header(page_number))
#v(2.5cm)
#weird(16pt, 11pt, "1. S", "UMMARY") 

The development of programs focused on chemistry has saved time, money, and, 
most importantly, reagents. Programs such as Gaussian, by solving Schrodinger's 
equation, have enabled a speed-up in chemical and technological development that 
was previously inconceivable. Other programs, such as LeaP, allow simulations of 
systems of hundreds of thousands of atoms, making it possible to simulate protein-drug 
interactions and obtain potential drugs for diseases that were previously incurable.

In order to solve the necessary equations, analytical solutions are not always 
possible, so so-called numerical methods are used. Numerical methods are mathematical 
techniques that, instead of finding an exact solution, as analytical methods do, 
find an approximate solution, which means they involve an error (the difference 
between the solution obtained and the exact solution) that is not necessarily 
negligible in the answer. The magnitude of the error depends on the numerical 
method, where the most accurate ones are usually the most expensive in terms of 
time and computation, so in some cases the method with the most error may be preferred.

All programs start from code in a file, and this TFG will explain the creation 
and development of a program in Python, where the program in question will attempt 
to solve equations posed in the subject “Physical Chemistry III” of the chemistry 
degree at the UB, focusing on the topic of electrochemistry. Relevant numerical 
methods for solving differential equations will also be mentioned, specifically the 
Crank-Nicolson method, an implicit method, and Runge-Kutta methods of order 2 and 4, 
which are explicit methods.

#v(3em)

*Keywords: * 

#pagebreak()

// NOTE: Resum

#let page_number = page_number + 1
#set page(header: which_header(page_number))
#v(2.5cm)
#weird(16pt, 11pt, "2. R", "ESUM") \
#v(1em)

El desarrollo de programas enfocados para la química ha permitido el ahorro de 
tiempo, dinero y, más importante, reactivos. Programas como gaussian, mediante 
la resolucion de ecuacion de Shrodinger, han permitido una agilicación en el 
desarrollo químico y tecnológico imposibles de concebir previamente.
Otros programas como LeaP, permiten simulaciones de sistemas de átomos del
orden de cientos de miles, pudiendo simular interacciones proteína-fármaco y 
obtener potenciales fármacos contra enfermedades antes sin cura.\ 
Para poder resolver las ecuaciones necesarias, las soluciones analíticas no 
siempre son posibles, sinó que se utilizan los llamados métodos numéricos. Los 
métodos númericos son técnicas matemáticas que en vez de encontrar una solución 
exacta, como lo hacen las analíticas, encuentras una solución aproximada, por lo 
que conllevan un error (la diferéncia entre la solución obtenida y la exacta) 
que no neceáriamente es negligible en la respuesta. La magnitud del error depende 
del método numérico, donde los más precisos suelen ser los más caros en cuanto a 
tiempo y computacionalmente, por lo que en algunos casos el método con más error 
puede ser el preferido. \

Todos los programas parten de código en un archivo, y en este TFG se explicará 
la creación y desarrollo de un programa en Python, donde el programa en cuestión 
tratará de resolver ecuaciones planteadas en la asignatura "Química Física III" 
del grado de química de la UB, centrado en el tema de electroquímica. También se 
mencionarán métodos numéricos relevantes en la resolución de ecuaciones 
diferenciales, específicamente el método de Crank-Nicolson, un método implícito, 
y Runge-Kutta de orden 2 y 4, métodos explícitos.

#v(3em)

*Palabras clave: * 

#pagebreak()

// NOTE: Introduction

#let page_number = page_number + 1
#set page(header: which_header(page_number))
#v(2.5cm)
#weird(16pt, 11pt, "3. I", "NTRODUCTION") 

 // WARN: Haciendose
== 3.1 Electrochemistry 
Electrochemistry is the branch of chemistry that focuses on reactions where 
oxidation and reduction happens. These reactions are commonly called redox 
reactions, inheriting the "red" from reduction and the "ox" from oxidation.
A general form for an oxidized species being reduced or vice-versa would be:
#v(0.5em)
$ O + e^-  harpoons.rtlb R $ \
#v(-1.5em)
Where $O$ represent the oxidized species, $R$ the reduced species, and $e^-$ 
stands for the electron particle with its negative charge. To express the 
reaction rate of a electrochemical reaction, the product of the reactive 
and a constant (cathodic or anodic) is taken, such that:
$ nu_c = k_c [O] \ nu_a = k_a [R] $ \
=== 3.1.1 The Butler-Volmer's equation
$ nu_c = k_c [O] \ nu_a = k_a [R] $ \
$ k_c [O] = A exp((Delta G^dagger.double_c)/(R T)) \ k_a [R] = A exp((Delta G^dagger.double_a)/(R T)) $
=== 3.1.2 Lievich and rotating electrode
=== 3.1.3 Battery discharge equation
=== 3.1.4 Cyclic Voltammetry
== 3.2 Numerical methods and algorithms used
=== 3.2.1 Euler's method
=== 3.2.1 RK2/Hund
=== 3.2.1 RK4
=== 3.2.1 Crank-Nicolson
=== 3.2.1 Thomas

#pagebreak()

// NOTE: Objectives

#let page_number = page_number + 1
#set page(header: which_header(page_number))
#v(2.5cm)
#weird(16pt, 11pt, "4. O", "BJECTIVES") 

// WARN: Objetivos muy simplistas

The primary objective of this TFG is pedagologycal and didactical: the creation, 
coding and implementation of the learnt electrochemical systems in class and 
the knowledge aquired through them is what's the true goal of this work.

#pagebreak()

// NOTE: Experimental section

#let page_number = page_number + 1
#set page(header: which_header(page_number))
#v(2.5cm)
#weird(16pt, 11pt, "5. E", "XPERIMENTAL SECTION") 

 // FIX: NO HECHO



#pagebreak()


// NOTE: First discussion title

#let page_number = page_number + 1
#set page(header: which_header(page_number))
#v(2.5cm)
#weird(16pt, 11pt, "6. F", "IRST DISCUSSION TITLE") 

 // FIX: NO HECHO



#pagebreak()

// NOTE: Second discussion title

#let page_number = page_number + 1
#set page(header: which_header(page_number))
#v(2.5cm)
#weird(16pt, 11pt, "7. S", "ECOND DISCUSSION TITLE") 

 // FIX: NO HECHO



#pagebreak()

// NOTE: Third discussion title

#let page_number = page_number + 1
#set page(header: which_header(page_number))
#v(2.5cm)
#weird(16pt, 11pt, "8. T", "HIRD DISCUSSION TITLE") 

 // FIX: NO HECHO



#pagebreak()

// NOTE: Conclusion

#let page_number = page_number + 1
#set page(header: which_header(page_number))
#v(2.5cm)
#weird(16pt, 11pt, "9. C", "onclusions") 

 // FIX: NO HECHO



#pagebreak()

// NOTE: References and notes

#let page_number = page_number + 1
#set page(header: which_header(page_number))
#v(2.5cm)
#weird(16pt, 11pt, "10. R", "ERENCES AND NOTES") 

 // FIX: NO HECHO



#pagebreak()

// NOTE: Acronyms

#let page_number = page_number + 1
#set page(header: which_header(page_number))
#v(2.5cm)
#weird(16pt, 11pt, "11. A", "CRONYMS") 

 // FIX: NO HECHO



#pagebreak()
