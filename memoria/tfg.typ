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
it's a manner of enhancing industrical technological capabilities. \
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

#weird(10pt, 8pt, "1. S", "UMMARY") \
#weird(10pt, 8pt, "2. R", "ESUMEN") \
#weird(10pt, 8pt, "3. I", "NTRODUCTION") \
#h(1.5em)3.1. Electrochemistry \
#h(3em)3.1.1. The butler-volmer's equation \
#h(3em)3.1.2. Levich and rotating disk electrode \
#h(3em)3.1.3. Battery discharge equation \
#h(3em)3.1.4. Cyclic voltammetry \
#h(1.5em)3.2. Numerical methods for differential equations \
#h(3em)3.2.1. Euler's method \
#h(3em)3.2.2. Heun's method \
#h(3em)3.2.3. Runge-kutta 4 \
#h(3em)3.2.4. Crank-nicolson \
#h(3em)3.2.5. Thomas algorithm \
#weird(10pt, 8pt, "4. O", "BJECTIVES") \
#weird(10pt, 8pt, "5. E", "XPERIMENTAL SECTION") \
#h(1.5em)5.1. Setup \
#h(1.5em)5.2. Interface \
#h(1.5em)5.3. Helper functions \
#h(1.5em)5.4. Databases \
#h(1.5em)5.5. Main functions \
#weird(10pt, 8pt, "6. I", "MPLEMENTATION OF THE CODE") \
#h(1.5em)6.1. ButlerVolmer.py \
#h(1.5em)6.2. RotatingDiskElectrode.py \
#h(1.5em)6.3. BatteryDischarge.py \
#h(1.5em)6.4. CyclicVoltammetry.py \
#weird(10pt, 8pt, "7. R", "ESULTS AND DISCUSSION") \
#weird(10pt, 8pt, "8. I", "MPROVEMENTS AND OPTIMIZATIONS") \
#weird(10pt, 8pt, "9.C", "ONCLUSIONS") \
#weird(10pt, 8pt, "10. R", "ERERENCES AND NOTES") \
#weird(10pt, 8pt, "11. A", "CRONYMS") \
#weird(10pt, 8pt, "A", "PPENDIX: PROGRAM CODE") \




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

*Keywords: * Mathematics, Chemical problems, Programming, Python.

#pagebreak()

// NOTE: Resum

#let page_number = page_number + 1
#set page(header: which_header(page_number))
#v(2.5cm)
#weird(16pt, 11pt, "2. R", "ESUMEN") \
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

*Palabras clave: * Matemáticas, Problemas Químicos, Programación, Python.

#pagebreak()

// NOTE: Introduction

#let page_number = page_number + 1
#set page(header: which_header(page_number))
#v(2.5cm)
#weird(16pt, 11pt, "3. I", "NTRODUCTION") 

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
where u is the concentration and at time $n$ and distance $i$. By defining 
$lambda = (D Delta t) / (2 Delta x^2)$, and separating by time step it can 
be rewritten as:
$ -lambda u^(n + 1)_(n + 1) + (1 + 2lambda) u^(n + 1)_i - lambda u^(n + 1)_i = 
  lambda u^n_(n + 1) + (1 - 2lambda) u^n_i + lambda u^n_i $
As the terms of the RHS are all known, the obtained row of the matrix has the 
following form:
$ mat(-lambda, 1 + 2 lambda, -lambda) vec(u^(n+1)_(i-1), u^(n+1)_i, u^(n+1)_(i+1)) =
vec(lambda u^n_(i-1), (1 - 2lambda) u^n_i, lambda u^n_(i+1)) $
This has the form $a x_(i-1) + b x_i + c x_(i + 1)  = d$. When extended to the whole 
matrix, the matrix obtained is a tridiagonal matrix, easily solved by the Thomas 
algorithm.

=== 3.2.5 Thomas algorithm
The Thomas algorithm replaces standard gaussian elimination by providing a 
more efficient manner of solving tridiagonal matrices, consisting of a forward sweep 
and back substitution.
$ mat(b_1, c_1, , , 0; a_2, b_2, c_2, ; , a_3, b_3, dots.down; 
, , dots.down, dots.down, c_(n-1); 0, , , a_n, b_n) vec(x_1, x_2, x_3, dots.v, x_n) = 
vec(d_1, d_2, d_3, dots.v, d_n) $
The forward sweep consists of eliminating every element of the $a$ diagonal by 
performing row elimination, using the previous row's $b$ as the pivot. By reducing 
each row into having one less unknown, the final row ends up with only one unknown, 
$x_n$ , which can be easily solved algebraically. Now that $x_n$'s value is known,
the $n-1_"th"$ row can be solved.This holds recursively, substituting the newly 
found value into the previous row until the first row is reached, solving the 
full vector of unknowns.

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

This TFG's program was made in python 3.13. The main structure of the program is 
shown by diagram 1:
#figure(
image("esquema.png"),
caption: [High-level structure overview]
)
== 5.1 Setup
A virtual enviroment is used, so generation of a python3 venv is required.
The generation of the virtual enviroment was done using the following shell 
command: 
$"python3 -m venv .venv && source .venv/bin/activate && pip -r requirements.txt" $
This shell command is actually three commands merged into one by the "and" operator
$"&&"$. First, an empty virtual enviroment is created, then the "source" command 
executes the .venv/bin/activate shell script and the context is shifted from the 
shell enviroment to the venv enviroment. Finally, pip, a python packet manager, 
installs the required packages listed in the requirements.txt file, which are:
numpy for efficient arrays, scipy for constants and matplotlib for plotting graphs.
Once the venv is created and has the required packages installed, the program can 
be executed by typing python3 interface.py or just ./interface.py.\
== 5.2 Interface
The python file interface.py is the entrypoint of the program, which loads and 
parses the databases into memory and displays the UI by printing text in the 
terminal. When executed, the program first displays a selection of options ranging 
from one to five, as shown in Figure 1, where each mode (excluding 5) executes 
one of the four main functions.\
#figure(
image("modes.png"),
caption: [program modes])
== 5.3 Helper functions
Helper functions can be thought of rogue functions, as they do not pertain to 
any specific section. Parsing functions, as "getElectrons" or "error", used 
to check if the input value is physically sensible are helper functions.\
== 5.4 Databases
Databases store large amounts of data in information. Two databases are used 
for this project as .csv files (comma separated values), one in StandardPotentials.csv, 
a database that stores half reactions and their corresponding standard reduction 
potential in volts, while the other, found in BatteryValues.csv relates the SOC of a 
Li-ion battery to its OCV. These databases are parsed at runtime, when the program 
starts, and are stored as python dictionaries.
== 5.5 Main functions
When selecting a mode, not including exit, the program will shift execution from 
interface.py to ButlerVolmer.py, RotatingDiskElectrode.py, BatteryDischarge.py 
or CyclicVoltammetry.py. All these files contain functions which ask for inputs 
and will plot a graph using the values given as inputs.

#pagebreak()


// NOTE: Butler-Volmer y RDE

#let page_number = page_number + 1
#set page(header: which_header(page_number))
#v(2.5cm)
#weird(16pt, 11pt, "6. I", "MPLEMENTATION OF THE CODE")

== 6.1 ButlerVolmer.py
The Butler Volmer equation, as described in Eq.10, some assumptions had to be 
made: First, the system is composed by a three electrodes, a working electrode, 
a reference electrode and a counter electrode. Second, only one ionic pair 
contributes to the output intensity and no successive reactions are contemplated. 
Third, to provide a meaningful graph, the function must be provided a range, 
which will be given by the overpotential. 

The Butler-Volmer function has the following structure:
#figure(
image("ButlerVolmerStructure.png", width: 50%),
caption: [structure overview of ButlerVolmer.py])

First, inputs are obtained through the error.py module, which consists of a 
huge helper function which takes as an input a string and returns either a 
float or a int. This function's structure consists of a match-case list 
(other languages might use switch-case instead), where if the input string 
finds a match, executes a try-catch which checks if the input is valid, and 
if it's not asks again until a valid input is provided. It should be mentioned 
that physical constants are not needed to be provided, as are assigned its 
proper value through the scipy.constants module. The final input, the half 
reaction, does undergo a special helper function, which checks if the provided 
half reaction exists in the StandardPotentials database and obtains the 
associated standard equilibrium potential. The following variables are 
obtiained as mentioned: temperature, concentration of both the reduced species 
and oxidized species, $j_0$, the applied potential and the reference electrode 
potential (using the Standard Hydrogen Electrode as 0V). A function named 
"getEeq" is called next, using two helper functions, "getElectrons" and 
"getStoichCoeffs", to calculate the equilibrium potential through the Nerst 
equation. The overpotential $eta$ is calculated and a numpy ndarray 
is created through numpy's linspace function, resulting in a ndarray from
$-abs(eta)$ to $abs(eta)$. Finally, using numpy's exp function, which 
applies a exponential to each element of an array without any need to loop, 
$j_a$ and $j_c$ are assigned, and the currents array obtained from $j_0 * (j_a - j_c)$
is used to plot a $E "vs" i$ graph with the "plotGraph" helper function and saved 
with whichever name the user wants.

== 6.2 RotatingDiskElectrode.py
In the same spirit as before, first some assumptions have to be explicited to 
determine the scope of the function. The function calculates the limiting 
current through the Levich equation, assuming the system uses a Rotating Disk
Electrode or RDE to allow a definite diffusion layer $delta$ to exist by providing 
laminar flow toward the electrode surfave. The reaction must be always diffusion
controlled, as the electron transfer should not be the limiting factor for the 
creation of the boundary layer. Finally, the solution uses a Newtonian fluid 
with constant viscosity. The structure of the function is the following:

#figure(
image("RotatingDiskElectrodeStructure.png", width: 50%),
caption: [structure overview of RotatingDiskElectrode.py])

Inputs are once again handled by the error helper function, ensuring that 
each physical parameter is valid before proceeding with the calculation, and
the following variables are obtained from the user: the number of electrons 
involved in the redox process, the diffusion coefficient, the bulk concentration, 
the viscosity of the electrolyte, the electrode area and the maximum rotation 
speed in RPM. 
A numpy ndarray for the rotation speed is created using numpy's linspace function, 
ranging from 100RPM to the input provided. The RPM array is converted into angular 
velocity $omega$ by the equivalence relation $omega = (2 pi)/60 "RPM"$. Then the
diffusion layer thickness $delta$ is calculated using Eq.13 and the limiting current 
$i_"lim"$ is computed for the entire range. Finally, the results are plotted as a
graph of $i_"lim"$ vs. RPM using the "plotGraph" helper function. 

#pagebreak()

// NOTE: Battery discharge


== 6.3 BatteryDischarge.py
The simulation of a battery discharge process, as formulated in the previous 
introduction section, relies on continuous tracking of the State of Charge 
(SOC) and its relation with the Open Circuit Voltage (OCV). The assumptions 
integrated into this functions are the following: First, the discharge occurs 
at a constant current $I$, with a low enough value that $R_"int" >> eta$ and 
the simplyfing Eq.16 to the following equation is possible: 
$ "SOC(t)" = "SOC"_(t=0) - I(t)/Q_"nom" $
Second, the internal resistance $R_"int"$ is treated as a constant parameter 
throughout the discharge process, as only one type of battery is considered, a 
Li-ion battery. Finally, the relation between the SOC and the OCV is derived 
from empirical data stored in the database BatteryValues.csv, interpolating 
linearly between known data points. The structure of the battery discharge 
simulation is the following:

#figure(
image("BatteryDischargeStructure.png", width: 50%),
caption: [structure overview of BatteryDischarge.py])

Inputs are once again handled by the error helper function. Variables as the 
current $I$, the internal resistance $R_"int"$, the nominal capacity $Q_"nom"$, 
initial SOC $"SOC"_(t = 0)$ and finally a time delimiter are obtained through the 
user's input. 
An empty python list is created for the output voltages, as unfortunately is 
impossible to take advantage of Numpy's ndarrays, as the end time is only known 
when the limiter is set: the dynamic duration of the battery's SOC is the source 
of the problem. Then, the SOC differential is calculated through the simplified 
Eq.37 ($d"SOC" = - I(t)/Q_"nom"$) and the if the time limiter is unset (input is 0) 
then it's set as 9e34 to ensure the battery discharges fully. The main loop of 
the function starts: for each step, the voltage through Eq.17 is 
calculated and appended into its corresponding position in the voltages list. 
The new SOC is set and the time differential is added to the current time, the 
loop is ready to perform a new iteration. To obtain the OCV atrributed to the 
current times' SOC, the "interpOCV" function is used. It's a function that takes 
the current SOC as an input and checks if the exact value exists in the database. 
If the value does not exist, it returns a value calculated through linear 
interpolation of the two closest values found in the database. As the condition of 
the while loop relies only on whether the time has reached or exceeded the delimited 
time or the SOC is a positive non-zero value, if the next time's SOC is negative no 
problem arises. Once the while loop loops through its last iteration, a graph is 
plotted through the "plotGraph" function using the voltages array turned into a 
Numpy ndarray with the array function and the time array, which is created using 
Numpy's "arange" with the length of the voltages multiplied by dt as an input.
After generating a voltage vs. time plot, the user is prompted for a filename, 
and the resulting graph is saved, providing a visual representation of the battery's 
discharge curve.

 // WARN: haciendose

== 6.4 CyclicVoltammetry.py
The simulation of a Cyclic Voltammetry (CV) is the most complex main funcion file in 
this TFG, as it uses three different numerical solvers for differential equations: 
Heun, RK4, and Crank-Nicolson. The "CyclicVoltammetry" is the first one that takes an 
input, the order of the RK method, if 0, Crank-Nicolson is used, if 2 or 4, RK2 (Heun)
or RK4 is selected instead. It simulates the concentration profiles of the oxidized 
species ($C_O$) and the reduced species ($C_R$) as the electrode potential sweeps 
linearly over time. To ensure physical accuracy, the following assumptions are made:
First, the solution is semi-infinite, meaning that the furthest concentrations do 
not sense the electrode. Second, the electron transfer follows Butler-Volmer kinetics 
at the electrode surface, and the spatial grid $x$ is large enough to prevent boundary 
effects. Third, the solution is homogenous at $t = 0$, having the same concentration 
at every point of the solution. Fourth and last, the support electrolyte does not 
interact, and all migration effects are negligible. The CyclicVoltammetry function 
is structured into several components in this fashion:

 #figure(
  image("CyclicVoltammetryStructure.png"),
  caption: [structure overview of CyclicVoltammetry.py])

As before, he "error" function is used to check for valid/sensible physical inputs.
The following variables are obtained from user input: starting potential, vertex 
potential, the scan rate $v$, number of cycles, the starting bulk concentration 
for both the oxidized $C_O"(bulk)"$and reduced species, the diffusion coefficient, 
the reaction's rate constant $k_0$, the cathodic symmetry factor. Through the "getEo"
helper function, the number of electrons are obtained, along the standard reduction 
potential of the input half reaction. The following variables are calculated based on 
the input obtained values: $t_"max"$, $x_"max"$, $d x$ and $d t$. The times array and 
potential signals array is generated using the "init_time_potential" function. This 
function creates a triangular potential waveform based on the input starting 
potential, vertex potential and scan rate ($v$). The function uses a phase shift and 
the modulo operator to ensure the potential array correctly ranges between the minimum 
potential $E_"start" - E_"vertex"$ and maximum potential $E_"start" + E_"vertex"$ for
the input number of cycles by the user. Then the space grid array $x$ is generated through
Numpy's "arange" function, always of size 301 values, where $x_0$ is the surface of the 
electrode and $x_301 = x_"max"$. Then the concentration profiles are built through 
Numpy's "full" function, filling two arrays of the same size than the times array with 
the input bulk concentrations for each species. Next, a python dictionary is created 
to provide all the needed variables to the functions in a less clutered manner, and 
finally, through if statements a numerical method is selected based on the order.

When order is 0, the Crank-Nicolson method is selected. The implementation of the 
the method employs a function named "newton_thomas" as a solver, which applies 
the implicit Crank-Nicolson method by averaging the current ($u_t$) and future 
($u_(t + 1)$ three-point formulas, as done in Eq. . Then solves the resulting 
tridiagonal system with the implemented Thomas algorithm. Because Crank-Nicolson 
cannot solve non-linear equations and the surface concentrations are non-linear, 
a Newton-Raphson iteration to find a self-consistent solution (a solution with 
low enough error) for $C_(O,0)$ and $C_(R,0)$ at each time step is employed.


When order is 2 is Heun's method (RK2) is selected. The implementation employs the 
a tecnique called Method of Lines, which discretizes the spacial second 
derivative of Eq. , effectively transforming the diffusion PDE into a system 
of coupled (ODEs). To solve the system obtained first, the algorithm calculates 
an initial derivative using the "get_derivatives" function and then computes an 
average slope over the interval $Delta t$ to as a slope corrector. This approach 
achieves a time precision of $O(Delta t^2)$.

Finally , when order is 4,the Runge-Kutta 4 (RK4) algorithm is employed. This 
method has the highest time precision of all the available methods in the file, 
obtaining derivatives at four distinct points within each time step: at the 
beginning ($k_1$), two estimates at the midpoint ($k_2, k_3$) and the end ($k_4$). 
This achieves a temporal truncation error of $O(Delta t^4)$. It is important to 
note that for all Runge-Kutta methods, the spatial precision remains $O(Delta x^2)$ 
due to the second-order central difference approximation being used for the diffusion 
term, which is the same for all methods of the Runge-Kutta family. 
Throughout the Heun and RK4 simulations, surface kinetics are handled by the 
"solve_surface_analytic" function. This function treats the electrode surface as  
a flux balance between the rate of diffusion and the Butler-Volmer electron transfer 
kinetics. 

Finally, for all the methods, the current density $j$ is calculated as the product
of the reaction rate and the Faraday constant. The concentration profiles for each
step are saved for animation purposees and the voltammogram is plotted as a 
standard $E$ vs. $j$ graph using the "plotGraph" helper function. The concentration 
profiles are also plotted as $C_O$ vs $x$ and $C_R$ vs $x$, using the profiles at 
time steps as frames for the animation.

#pagebreak()

#let page_number = page_number + 1
#set page(header: which_header(page_number))
#v(2.5cm)
#weird(16pt, 11pt, "7. R", "ESULTS AND DISCUSSION") 

To prove the program works and the simulations are legitimate, the graphical output 
from the simulations must be compared with experimental sources. 

== 7.1 Butler Volmer
The Butler volme output graph must imitate closely the 

#pagebreak()

#let page_number = page_number + 1
#set page(header: which_header(page_number))
#v(2.5cm)
#weird(16pt, 11pt, "8. I", "MPROVEMENTS AND OPTIMIZATIONS") 

The current program can be improved and optimized further, however, as the program 
was made for a pedagological point some limits were assumed. This section will 
cover further improvements for the overall structure and performance of the program.\
First, the UI is very minimal, only a few lines of text printed to the terminal. 
The UI could be improved greatly by using tools like Streamlit or Custom Tkinter.
Streamlit could be used to turn the program into a web app, providing it of a simple 
yet elegant interface, while Custom Tkinter would make of the program an application, 
with its own window and a personalized interface, tailored to the code written.\
To make the code more readable, Object Oriented Programming (OOP), along classes and 
methods could make the code modular and more readable, in case other programmers 
join the codebase. Other external modules would too make the code more compact, as 
other than pedagologically, no real need to write the numerical algorithms existed. 
Using already existing functions from well known modules will make the code clearer 
and possibly faster, as several optimizations can be done for the solvers. \ Another 
drawback is the way inputs are provided to the program. Thanks to the "error" function 
which handles bad inputs the program won't crash, but for batches of simulations, 
manually inputing each parameter every time a simulation is needed might be a cause of 
unnecessary burden for the user. Taking as input a .json file or a .yaml would be a 
gret improvement over the current input system. Although less intuitive for new users, 
the easy reproductibility and the possibility of providing multiple inputs would allow 
for easy parallelization and for supercomputers to not wait idle waiting while another 
input is being sent.\
On the topic of parallelization, python 3.14 removes the Global Interpreter Lock (GIL)
allowing for true multithreading in Python 3. This update allows for the python
interpreter not to be locked into only being used by one thread. By rewritting the 
program for Python 3.14, and importing the new multithreading module, making correct 
use of threading, an exponential performance improvement on current code would be 
observed. Without any version change, the program can be made concurrent by Inter 
Process Comunication (IPC), creating new processes and comunicating them. Although 
more expensive than multithreading, can be achieved without any version change. \
Python's interpreter is a massive slowdown compared to compiled languages, rewritting 
the project in other programming languages which are compiled, like C, C++, or new 
programming languages like Zig or Rust would benefit from compiler optimizations, 
enabling a huge optimization. It should be mentioned that if performance is the 
only objective even if compiled, programming languages like Java or C\# should be 
avoided, as not requiring memmory mangement implies a underlying garbage collector, 
which slows down the runtime process. To truly optimize it fully, the cache should 
be taken advantage of as much as possible. Compressing the data needed in such a 
way all is contained in a the cache will turn the program into a very performant 
simulator. Avoiding cache flushes can lessen the runtime by orders of magnitude, 
especially if memory is stored in L1 cache or L2 instead. Albeit buying CPUs with 
megabytes of cache can improve cache hits, optimizing the code will be the principal 
manner to achieve true performance. \ Finally, vectorization through Single Instruction, 
Multiple Data (SIMD) coupled with branchless programming will take the code to another 
level of performance, as no branch predictions will be needed for the CPUs and SIMD 
instructions perform multiple operations in one CPU cycle. 


#pagebreak()

// NOTE: Conclusion

#let page_number = page_number + 1
#set page(header: which_header(page_number))
#v(2.5cm)
#weird(16pt, 11pt, "9. C", "ONCLUSIONS") 

To conclude this TFG, the author would like to think the pedagological objective has 
been achieved. A program has been made, which simulates accurately real phenomena 
using advanced numerical methods, as seen in the sections before.
#pagebreak()

// NOTE: References and notes

#let page_number = page_number + 1
#set page(header: which_header(page_number))
#v(2.5cm)
#weird(16pt, 11pt, "10. R", "ERENCES AND NOTES") 

1. Bard, A. J., & Faulkner, L. R. (2001). Electrochemical Methods: Fundamentals and Applications. John Wiley & Sons.\
2. Atkins, P., & de Paula, J. (2014). Physical chemistry (10th ed.). Oxford University Press.
3. Rumble, J. R. (Ed.). (2021). CRC handbook of chemistry and physics (102nd ed.). CRC Press.
4. Zhang, S., Zhang, Y., Liu, Y., & Zhang, J. (2022). SOC estimation of lithium-ion battery based on an improved equivalent circuit model and an improved adaptive unscented Kalman filter. Energies, 15(23), 9142. https://doi.org/10.3390/en15239142
5. Press, W. H., Teukolsky, S. A., Vetterling, W. T., & Flannery, B. P. (2007). Numerical recipes: The art of scientific computing (3.ª ed.). Cambridge University Press.
6. Elgrishi, N., Rountree, K. J., McCarthy, B. D., Rountree, E. S., Eisenhart, T. T., y Dempsey, J. L. (2018). A Practical Beginner’s Guide to Cyclic Voltammetry. Journal of Chemical Education, 95(2), 197–206 . https://doi.org/10.1021/acs.jchemed.7b00361.

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
- UI: User Interface
- LHS: Left Hand Side
- RHS: Right Hand Side
- parse: analyse text data and transform it to manageable data 

#pagebreak()

#let page_number = page_number + 1
#set page(header: which_header(page_number))
#align(center)[#h(5cm)#set text(fill: black.lighten(35%))
#weird(36pt, 26pt, "A", "PPENDICES")]

#pagebreak()

#let page_number = page_number + 1
#set page(header: which_header(page_number))
#weird(16pt, 11pt, "P", "ROGRAM CODE")
#let page_number = page_number + 1
#set page(header: which_header(page_number))
#image("interface.png")
#let page_number = page_number + 1
#set page(header: which_header(page_number))
#image("utils.png")
#let page_number = page_number + 1
#set page(header: which_header(page_number))
#image("error1.png")
#let page_number = page_number + 1
#set page(header: which_header(page_number))
#image("error2.png")
#let page_number = page_number + 1
#set page(header: which_header(page_number))
#image("error3.png")
#let page_number = page_number + 1
#set page(header: which_header(page_number))
#image("ButlerVolmer.png")
#let page_number = page_number + 1
#set page(header: which_header(page_number))
#image("RotatingDiskElectrode.png")
#let page_number = page_number + 1
#set page(header: which_header(page_number))
#image("BatteryDischarge.png")
#let page_number = page_number + 1
#set page(header: which_header(page_number))
#image("CyclicVoltammetry1.png")
#let page_number = page_number + 1
#set page(header: which_header(page_number))
#image("CyclicVoltammetry2.png")
#let page_number = page_number + 1
#set page(header: which_header(page_number))
#image("CyclicVoltammetry3.png")
#let page_number = page_number + 1
#set page(header: which_header(page_number))
#image("CyclicVoltammetry4.png")


/* INFO:        NOMBRE           |       HECHO
            agradecimientos               V
                  SDG                     V
                Summary                   V
                 Resum                    V
              Introduccion         V V V V V V V V 
              Estructura                  V
               Ex 1 y 2                  V V
                 Pilas                    V
                  CV                      V
               Conclusion               V X v
*/
