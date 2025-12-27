#set document(title: [Manica_Georigiev_TFG])
#set math.equation(numbering: "[Eq.1]", number-align: left)
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
most importantly, reagents. Programs such as Gaussian, by solving Schrödinger's 
equation, have enabled a speed-up in chemical and technological development that 
was previously inconceivable. Other programs, such as LeaP, allow simulations of 
systems of hundreds of thousands of atoms, making it possible to simulate protein-drug 
interactions and obtain potential drugs for diseases that were previously incurable.

In order to solve the necessary equations, analytical solutions are not always 
possible, so so-called numerical methods are used. Numerical methods are mathematical 
techniques that, instead of finding an exact solution as analytical methods do, 
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
la resolucion de ecuacion de Schrödinger, han permitido una agilicación en el 
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
$ O + e^-  harpoons.rtlb R $ 
A cathode is an electrode where the reduction process ocurrs, while an anode 
is an electrode where an oxidation process happens. Electrochemical reactions 
and processes might be referred as anodic or cathodic instead. To express the 
reaction rate of a electrochemical reaction, the product of the reactive and 
a constant (cathodic or anodic) is taken, such that:
$ nu_c = k_c [O] \ 
  nu_a = k_a [R] $ 
The constants, $k_a$ and $k_c$, can be represented through the transition state 
theory as:
$ k_c = A exp((Delta G^dagger.double_c)/(R T)) \ 
  k_a = A exp((Delta G^dagger.double_a)/(R T)) $
Before applying any potential, $k_a = k_c$ since $Delta G^dagger.double_a = Delta 
G^dagger.double_c$, so $k_0$ is defined to be either $k_a$ or $k_c$ at $E = E_0$.
When applying a potential to an electrode, the energy of the electrons in it canges,
which shifts the activation energy barriers. The shift can be represented by the 
transfer coefficient ($beta$), which represents the fraction of the electrical energy 
that affects the transition state. Since $E^0$ is the formal potential, the change in 
activation energy due to the overpotential relating to it ($eta = E_(a p) - E^0$) is:
$ Delta G^dagger.double_c (E) = Delta G^dagger.double_(0,c) + beta F (E - E^0) \
 Delta G^dagger.double_a (E) = Delta G^dagger.double_(0,a) + (1 - beta) F (E - E^0) $
The cathodic barrier shifts into higher energies when the overpotential is positive 
by a factor of $beta F eta$, while the anodic shifts into lower energies by a factor 
of $(1 - beta) F eta$. The opposite is true for negative overpotentials. Substituting 
brings this new form of Eq.3 alive:
$ k_c = k_0 exp((Delta G^dagger.double_c)/(R T)) \ 
  k_a = k_0 exp((Delta G^dagger.double_a)/(R T)) $

=== 3.1.1 The Butler-Volmer's equation
The Butler-Volmer equation connects the thermodynamics of the reactions, given by 
the overpotential ($eta$), which acts as the driving force of the reaction, with 
kinetics given by the total current $i$, which is a sum of the partial anodic and 
cathodic currents $i_a$ and $i_c$. 
The common writing of the Butler-Volmer uses current densities, which are defined as:
$ j = i / A $
Since $j = j_a - j_c$ and $j = F nu$, then follows $j = F (nu_a - nu_c)$ and using 
Eq.2, an expanded version can be written as:
$ j = F k_0 ([R]_s k_a - [O]_s k_c) $
The resulting equation uses the concentrations on the surface of the electrode, as 
there's where the reaction happens, and it can be further expanded using Eq.3:
$ j = F k_0 ([R]_s exp(((1 - beta) F eta)/(R T)) - [O]_s exp((beta F eta)/(R T))) $
The final expanded equation now relates the overpotential directly with the current 
density. This, however, is not the common writing of the Butler-Volmer equation 
either, as it can be compacted by defining the exchange current density $j_0$ as:
$ j_0 = F k_0 [R]^beta [O]_b^((1 - beta)) $
$j_0$ is defined as he dynamic equilibrium current that flows equally in both 
directions when the overpotential is zero, ie. $eta = 0$, and reflects the 
intrinsic speed of the reaction. The exchange current density uses the bulk 
concentrations as no reaction is happening, so the system is at rest. Compacted, 
the Butler-Volmer equation is written commonly as: 
$ j = j_0 (exp(((1 - beta) n F eta)/(R T)) - exp(-(beta n F eta)/(R T))) $

=== 3.1.2 Levich and rotating disk electrode
A Rotating Disk Electrode (RDE) is a tool used to measure the hydrodynamics and 
study reaction kinetics and mass transport. Unlike a stationary electrode, 
diffusion does not depend on a stagnant layer that grows over time, since through 
rotation the RDE creates a constant, well-defined flow of electrolyte toward the 
surface of the electrode.\
The Levich equation is able to link the rotation of the RDE with the resulting 
limit current $i_L$, which is defined as the maximum current achieved when the 
reaction is diffusion controlled, ie., by how fast reactants can be diffused 
through the boundary layer, throughout this formula:
$ i_L = 0.620n F A D^(2/3) omega^(1/2) nu^(-1/6) C $
which can be written alternatively, integrating the diffusion layer as a parameter: 
$ i_L = n F A D C / delta $
with the diffusion layer being defined as:
$ delta = 1.61 D^(1/3) omega^(-1/2) nu^(1/6) $

=== 3.1.3 Battery discharge equation
Voltage though a circuit can be measured and calculated with Ohm's law:
$ V = I R $
A battery can be thought as a simple circuit, with an internal resistance:
$ V = V - I R_"int" $
The voltage when there's no current is known as the open circuit voltage, and 
is commonly denoted by the initials OCV, and depends on the state of charge of 
the battery (SOC), which is a percentage that represents the capacity in a 
rechargeable battery, where at 100% the battery is at full capacity and at 0% 
the battery is depleted. The discharge of a battery can be thought as the 
decreasing of the SOC. The SOC, then, can be thought of as a function depending 
on time, where the current $I(t)$ acts as a drain to the capacity ($Q_"nom"$) of 
the battery, expresed as:
$ "SOC"(t) = "SOC"_(t = 0) - 1 / Q_"nom"integral^t_0 I(t) d t $
The voltage of the battery then can be expressed as the difference of the OCV 
and the sum of causes of energy loss:
$ V = "OCV"("SOC") - I R_"int" - eta_"act" $
While $I R_"int"$ accounts for energy loss as a "physical" resistance, as electrons 
are moving through the material, $eta_"act"$ accounts for energy loss as a "chemical" 
resistance, as the energy required to make the redox reaction happen. Additionally, 
$eta_"act"$ obeys the following relation:
$ eta_"act" =  (R T)/(beta F) ln (I / I_0) $
Thus, for low currents $I R_"int" >> eta_"act"$, and the voltage of a battery can be 
simplified as the following equation:
$ V = "OCV"("SOC") - I R_"int" $

=== 3.1.4 Cyclic Voltammetry
A Cyclic Voltammetry is a electrochemical measuring technique where electrical 
potential $E$ is applied in sweeps and current is obtained as an output. The current 
$i$ obtained is determined by the Butler-Volmer equation (Eq.10), however, another 
important factor to consider is the diffusion of the electrolyte, which is dictated 
by Fick's laws:
$ J = - D (partial C) / (partial x) $
The first fick law describes how a diffusion flux is formed when a concentration 
gradient exists throughout a distance.The negative sign is crucial because it 
shows that the species move from high to low density regions to reach equilibrium.
$ (partial C) / (partial t) = D (partial^2 C) / (partial x^2) $
The second fick law is obtained by differentiating the first law and describes 
how the rate of change of concentration respect to time relates to the curvature 
of the concentration gradient.\
It is useful to think of the CV sweep as a transition between two regimes. 
First, the system is kinetically controlled, as the Butler-Volmer equation 
dictates the current based on the applied overpotential. However, as no infinite  
concentration of electrolyte exists, when the overpotential exceeds a threshold, 
the concentration at the surface of the electrode drops to zero and the system 
becomes diffusion controlled. At that point, diffusion of the electrolyte to the 
surface of the electrode becomes the limiting factor and the current begins to 
decrease as the diffusion layer extends further into the bulk solution.

== 3.2 Numerical methods for differential equations
Numerical methods for differential equations are specific numerical methods for 
solving differential equations, which are equations that relate an unknown function 
to its derivatives. Some numerical methods approximate derivatives with finite 
differences. Finding a first order derivative is as simple as substituting the 
infinitesimally small $h$ by a finitely small $Delta x$ in the definition of 
derivative:
$ f'(x) approx (f(x + Delta x) - f(x)) / (Delta x) $
This substitution cannot obtain the same value as the real derivative, and 
brings a truncation error of the order of $Delta x$, commonly denoted $O(Delta x)$. 
For higher order derivatives, we can repeat the same process, however, the 
associated error can be minimized by using a clever trick with Taylor expansions:
$ f(x + Delta x) = f(x) + f'(x) Delta x f''(x) Delta x^2 + f'''(x) Delta x^3 + O(Delta x^4) $
$ f(x - Delta x) = f(x) - f'(x) Delta x f''(x) Delta x^2 - f'''(x) Delta x^3 + O(Delta x^4) $
Adding $f(x + Delta x) + f(x - Delta x)$ we obtain the three point formula:
$ f''(x) approx (f(x + Delta x) - 2f(x) + f(x - Delta x)) / (Delta x^2) $
This numerical second order derivative has a error of order $O(Delta x^4)$, as 
both the first and third order terms vanish when adding, which makes it both 
efficient and precise for numerical computations.

=== 3.2.1 Euler's method
The simplest numerical method to solve differential equations was published by 
Leonhard Euler in 1768, and it's known as Forward Euler, as it's an explicit
method, or just Euler's method. This numerical method core idea is to treata 
differential equation as a formula for the slope of the tangent line. The curve 
$f(x)$ is initially unknown, but a defined starting point $(x_0, y_0)$ is required,
as computing the derivative $f'(x)$ at this point enables the recursive calculation 
of the next points. The tangent line obtained through the derivative acts as the 
direction towards the next step, the distance moved depending on $Delta x$:
$ f(x_0 + Delta x) = f(x_0) + f'(x_0) Delta x $
or written alternatively, using subscripts instead:
$ y_1 = y_0 + y'_0 Delta x \
  y_2 = y_1 + y'_1 Delta x \
... $
Recursively, we arrive to a general formula:
$ y_(n + 1) = y_n + y'_n Delta x $

=== 3.2.2 Heun's method
Heun's method, also named improved forward Euler's method fixes Euler's method 
by enforcing a slope correction, it being the average of the current slope and the 
Euler's predicteed next point. It adresses the primary weakness of Euler's approach, 
the assumption that the slope remains constant over the entire interval $Delta x$.
The process begins by calculating an initial estimate of the next point using a 
standard Euler step. Although more expensive to compute, averaging this new slope 
with the initial slope brings a much more accurate result, as it accounts for the 
curvature of the function $f(x)$. The resulting function is written as:
$ y_(n+1) = y_n + (Delta x)/2 (k_1 + k_2) $
where the slopes $k_1$ and $k_2$ are:
$ cases(k_1 = f(x_n, y_n),
  k_2 = f(x_n + Delta x, y_n + Delta x k_1)) $
where $f(x_n, y_n)$ corresponds to the derivative at point $x_n$ given the value 
$y_n$.

=== 3.2.3 Runge-Kutta 4
The Runge-Kutta methods are a family of explicit numerical methods, as it's a 
generalized form for explicit ODE methods. Euler's method is taken as base, and 
slope corrections are introduced, with higher order Runge-Kutta methods being more 
precise, at the cost of being more computationally expensive. Euler's method by 
definition is RK1, then Heun's method is RK2, with the most famous method being 
the fourth order RK4. RK4 expands on Heun's slope correction: 
$ y_(n+1) = y_n + (Delta x)/6 (k_1 + 2k_2 + 2k_3 + k_4) $
with slopes corresponding to:
$ cases(
  k_1 = f(x_n, y_n),
  k_2 = f(x_n + (Delta x)/2, y_n + (Delta x)/2 k_1),
  k_3 = f(x_n + (Delta x)/2, y_n + (Delta x)/2 k_2),
  k_4 = f(x_n + Delta x, y_n + Delta x k_3)) $

=== 3.2.4 Crank-Nicolson
The Runge-Kutta family of methods provide high accuracy for ordinary differential 
equations, however, for partial differential equations, which involve more than 
one variable, significant stability challenges are introduced.
The Crank-Nicolson method is a second-order implicit method designed to overcome 
the numerical instability that arises with explicit methods when the time step 
$Delta t$ is large relative to the distance step $Delta x$. It uses a finite 
difference method for the derivatives, using the three point formula, to solve 
differental equations similar to Fick's second law, where a first order derivative 
in the LHS with respect to one variable is related to a second order derivative in 
the RHS with respect to another variable. Then it averages the current time step 
with the next time step.
Using the Crank-Nicolson for the diffusion equation (Eq.21) produces the following 
equation:
$ (u^(n + 1)_i - u^n_i) / (Delta t) = D/2 ((u_(i + 1)^n - 2u_i^n + u_(i - 1)^n) / (Delta x^2) +  
(u_(i + 1)^(n + 1) - 2u_i^(n + 1) + u_(i - 1)^(n + 1)) / (Delta x^2)) $
where u is the concentration and at time n and distance i. By defining 
$r = (D Delta t) / (2 Delta x^2)$, and separating by time step it can 
be rewritten as:
$ -r u^(n + 1)_(n + 1) + (1 + 2r) u^(n + 1)_i - r u^(n + 1)_i = 
  r u^n_(n + 1) + (1 - 2r) u^n_i + r u^n_i $
As the terms of the RHS are all known, the obtained row of the matrix has the 
following form:
$ mat(-r, 1 + 2 r, -r) vec(u^(n+1)_(i-1), u^(n+1)_i, u^(n+1)_(i+1)) =
vec(r u^n_(i-1), (1 - 2r) u^n_i, r u^n_(i+1)) $
This has the form $a x_(i-1) + b x_i + c x_(i + 1)  = d$. When extended to the whole 
matrix, the matrix obtained is a tridiagonal matrix, easily solved by the Thomas 
algorithm.

=== 3.2.5 Thomas algorithm
The Thomas algorithm replaces standard gaussian elimination by providing a 
more efficient way of solving tridiagonal matrices, consisting of a forward sweep 
and back substitution.
$ mat(b_1, c_1, , , 0; a_2, b_2, c_2, ; , a_3, b_3, dots.down; 
, , dots.down, dots.down, c_(n-1); 0, , , a_n, b_n) vec(x_1, x_2, x_3, dots.v, x_n) = 
vec(d_1, d_2, d_3, dots.v, d_n) $
The forward sweep consists of eliminating every element of the $a$ diagonal by performing 
row elimination, using the previous row's $b$ as the pivot. By reducing each row into 
having one less unknown, the final row ends up with only one unknown, $x_n$ , which can be 
easily solved algebraically. Now that $x_n$ has a value, the $n-1_"th"$ row can be solved.
This holds recursively, substituting the newly found value into the previous row until the 
first row is reached, solving the full vector of unknowns.

#pagebreak()

// NOTE: Objectives

#let page_number = page_number + 1
#set page(header: which_header(page_number))
#v(2.5cm)
#weird(16pt, 11pt, "4. O", "BJECTIVES") 

// WARN: Objetivos muy simplistas

The primary objective of this TFG is both pedagologycal and didactical: the 
creation, coding and implementation of the learnt electrochemical systems 
in class and the knowledge aquired through them are the true goal of this work.

#pagebreak()

// NOTE: Experimental section

#let page_number = page_number + 1
#set page(header: which_header(page_number))
#v(2.5cm)
#weird(16pt, 11pt, "5. E", "XPERIMENTAL SECTION") 

The structure of t


#pagebreak()


// NOTE: Butler-Volmer y RDE

#let page_number = page_number + 1
#set page(header: which_header(page_number))
#v(2.5cm)
#weird(16pt, 11pt, "6. B", "UTLER-VOLMER AND RDE") 

 // FIX: NO HECHO



#pagebreak()

// NOTE: Pilas de litio

#let page_number = page_number + 1
#set page(header: which_header(page_number))
#v(2.5cm)
#weird(16pt, 11pt, "7. B", "ATTERY DISCHARGE") 

 // FIX: NO HECHO



#pagebreak()

// NOTE: Voltamperometria ciclica

#let page_number = page_number + 1
#set page(header: which_header(page_number))
#v(2.5cm)
#weird(16pt, 11pt, "8. C", "YCLIC VOLTAMMETRY") 

 // FIX: NO HECHO



#pagebreak()

// NOTE: Conclusion

#let page_number = page_number + 1
#set page(header: which_header(page_number))
#v(2.5cm)
#weird(16pt, 11pt, "9. C", "ONCLUSIONS") 

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

- $[O]:$ oxidized species
- $[R]:$ reduced species
- $e^-:$ electron particle with its negative charge 
- $nu_a:$ anodic reaction rate 
- $nu_c:$ cathodic reaction rate
- $k_a:$ anodic constant 
- $k_c:$ cathodic constant 
- $A:$ Arhenius preexponent
- $Delta G^dagger.double:$  Gibbs activation energy
- $R:$ ideal gas constant
- $T:$ temperature
- $j:$ current density
- $j_0:$ exchange current density
- $i:$ current
- $I:$ current
- $I_0:$ exchange current
- $E_0:$ formal potential
- $n:$ number of electrons of the reaction
- $A:$ area
- $D:$ diffusion coefficient
- $J:$ diffusion flux
- $O(Delta x^n):$ Error of order n
- RK: Runge-Kutta
- ODE: Ordinary Differential Equation
- PDE: Partial Differential Equation

 // WARN: Haciendose

#pagebreak()


/* INFO:        NOMBRE           |       HECHO
            agradecimientos               V
                  SDG                     V
                Summary                   V
                 Resum                    V
              Introduccion         V V V V V X X X 
              Estructura                  X
               Ex 1 y 2                  X X
                 Pilas                    X
                  CV                      X
               Conclusion               X X X
*/
