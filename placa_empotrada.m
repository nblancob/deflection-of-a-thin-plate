function placa2 (n,m)
% Función para encontrar por diferencias finitas la deflexión de una placa
% delgada cuya ecuación diferencial esta dada por: (∂^4 w)/(∂x^4 )+2 (∂^4
% w)/(∂x^2*∂y^2 )+(∂^4 w)/(∂y^4 )=(q_z (x,y))/D, con parametros n=nodos en
% dirección x y m=nodos en dirección y.
%---------------------------------------------------------------------------------------------------------------------  
% Condiciones de la placa 
  q= 10000; % Carga distribuida en N/m^2
  E= 200e9; % Modulo de elasticidad en Pascales
  v= 0.3; % Relación de Poisson
  h=0.1; % Espesor de la placa en m
  D= (E*h^3)/(12*(1-v^2)); %Constante de rigidez a la flexión de la placa
%----------------------------------------------------------------------------------------------------------------------
% Condiciones de la malla
  X = 4; % Ancho de la placa
  Y = 3; % Largo de la placa
  hx = X/(n-1); % Espaciamiento en dirección x
  hy = Y/(m-1); % Espaciamiento en dirección y
  nnin = n-2; % Numero de nodos internos en dirección x
  nmin = m-2; % Numero de nodos internos en dirección y
 %----------------------------------------------------------------------------------------------------------------------
 % Condiciones de borde para  la placa dadas por la funcion borde que dependeran del tipo de apoyo, y
 % guardadas en la matriz c
 for i=1:n
    x = (i-1)*hx;
    c(i,1) = borde(x, 0, X, Y); % Para los condiciones de borde simplemente apoyadas o empotradas el valor de los extremos sera u=0 y z=0
    c(i,m) = borde(x, Y, X, Y);% En este caso los cuatro bordes se encuentran simplemente apoyados
  end 
  for j = 2:m-1
    y = (j-1)*hy; 
    c(1,j) = borde(0, y, X, Y);
    c(n,j) = borde(X, y, X, Y);
  end 
  %----------------------------------------------------------------------------------------------------------------------
  % Construcción de la matriz de coeficientes dados por la ecuación en
  % diferencias finitas:
  % u_(i+1,j)+u_(i-1,j)+u_(i,j+1)+u_(i,j-1)-〖4u〗_(i,j)=(h^2 q_z (x,y))/D
  % note que esta expresión esta en función de u y no de w, por lo que se
  % tendra que resolver un segundo sistema para encontrar en valor de la
  % deflexión
  kb = zeros(nnin); %Matriz de nodos internos
  kb(1,1) = -4;
  for i=1:n-3
    kb(i,i+1)=1 ;
    kb(i+1,i)=1 ;
    kb(i+1,i+1)=-4 ; 
  end 
  I = eye(n-2) ; 

  k = zeros(nnin*nmin) ; 
  k(1:nnin , 1:nnin) = kb ;
  for j=2:nmin 
    jj = nnin*(j-1)+1;
    k(jj:jj+nnin-1 , jj:jj+nnin-1) = kb;
    k(jj-nnin:jj-1 , jj:jj+nnin-1) = I;
    k(jj:jj+nnin-1 , jj-nnin:jj-1) = I;
  end 

  %-----------------------------------------------------------------------------------------------------------------------
  % Vector de terminos independientes dado por las condiciones de borde y
  % la expresión: (h^4 q_z (x,y))/D, donde h es igual al espaciamiento de
  % los nodos
  f = zeros(nnin*nmin , 1);
  for i = 1:nnin*nmin
    f(i) = (q*hx^2/D); 
  end 
%---------------------------------------------------------------------------------------------------------------------------
% Solución del sistema matricial obtenemos los valores para u
  u = k\f;
  u=u*hx^2;
%---------------------------------------------------------------------------------------------------------------------------
% Solución del segundo matricial para encontrar el valor de la deflexión w
% con la ecuación diferencial u=(∂^2 w)/(∂x^2 )+(∂^2 w)/(∂y^2 ) que
% puede ser expresada en diferencias como:  
   w=k\-u;
%---------------------------------------------------------------------------------------------------------------------------
% Inclusión de los bordes en la solución
  inud = 1;
  for i=2:n-1
    for j=2:m-1
      c(i,j) = w(inud) ;
      inud = inud + 1 ;
    end 
  end 
  
%----------------------------------------------------------------------------------------------------------------------------
% Representación grafica de la placa deflectada.
  x = 0:hx:X ;
  y = 0:hy:Y ;
  colormap("winter");
  surfl(x, y, c'),title('Grafica deflexión placa')%,shading("interp")
  a=min(c);
  dmax=min(a)  %Deflexión maxima
  end 