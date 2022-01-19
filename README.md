# deflection-of-a-thin-plate
En la ingeniería estructural se considera como placa al sólido tridimensional el cual tiene una de sus dimensiones, el espesor, mucho menor que las otras, el cálculo de placas tiene su origen con los trabajos que realizo Euler en el siglo XIII, a partir de los cuales se han desarrollado la teoría fundamental de placas. 
Para este trabajo se analizará los esfuerzos en una placa delgada mediante la teoría clásica de placas de Kirchhoff la cual desprecia la deformación por córtate, lo que permite expresar la ecuación con derivadas parciales en función de la deflexión.
 En esta ecuación se plantea la deflexión transversal w(x,y)  como una ecuación de cuarto orden, la cual requiere satisfacer dos condiciones de borde en cada extremo:


∇^4 w(x,y)=(q_z (x,y))/D      |      D=(Eh^3)/(12 (1-v^2))  
      
Para resolver esta ecuación existen diversos métodos aproximados en este caso se buscará solucionar la ecuación por el método de diferencias finitas implementado en el lenguaje de programación Matlab y explicado en el documento pdf.
