function placa (n,m)
% Función para encontrar por diferencias finitas la deflexión de una placa
% delgada cuya ecuación diferencial esta dada por: (∂^4 w)/(∂x^4 )+2 (∂^4
% w)/(∂x^2*∂y^2 )+(∂^4 w)/(∂y^4 )=(q_z (x,y))/D, con parametros n=nodos en
% dirección x y m=nodos en dirección y.
%---------------------------------------------------------------------------------------------------------------------  
% Condiciones de la placa 
  q=10000; % Carga distribuida en N/m^2
  E=200e9; % Modulo de elasticidad en Pascales
  v= 0.3; % Relación de Poisson
  h=0.1; % Espesor de la placa en m
  D= E*h^3/(12*(1-v^2)); %Constante de rigidez a la flexión de la placa
%----------------------------------------------------------------------------------------------------------------------
% Condiciones de la malla
  X = 4; % Ancho de la placa
  Y = 3; % Largo de la placa
  hx = X/(n-1);  % Espaciamiento en dirección x
  hy = Y/(m-1); % Espaciamiento en dirección y
  nnin = n-2; % Numero de nodos internos en dirección x
  nmin = m-2; % Numero de nodos internos en dirección y
 %----------------------------------------------------------------------------------------------------------------------
 % Condiciones de borde para  la placa dadas por la funcion borde que dependeran del tipo de apoyo, y
 % guardadas en la matriz c
 for i=1:n
    x = (i-1)*hx;
    c(i,1) = borde(x, 0, X, Y);
    c(i,m) = borde(x, Y, X, Y);
  end 
  for j = 2:m-1
    y = (j-1)*hy; 
    c(1,j) = borde(0, y, X, Y);
    c(n,j) = borde(X, y, X, Y);
  end 
  %----------------------------------------------------------------------------------------------------------------------
  % Construcción de la matriz de coeficientes dados por la ecuación en
  % diferencias finitas: 20w_(i,j)-8(w_(i+1.j)+w_(i-1,j)+w_(i,j+1)+w_(i,j-1) )+2(w_(i+1,j+1)+w_(i-1,j+1)+w_(i+1,j-1)+w_(i-1,j-1) )+w_(i+2,j)+w_(i-2,j)+w_(i,j+2)+w_(i,j-2)=(h^4 q_z (x,y))/D
  kb = ones(nnin); %Matriz de nodos internos
  kb(1,1) = 20;
  for i=1:n-3
    kb(i,i+1)=-8 ;
    kb(i+1,i)=-8 ;
    kb(i+1,i+1)=20; 
  end 
  kb(fix(n/2),fix(n/2))=19;
  kc=kb;
  kc(1,1)=21;
  for i=1:n-3
  kc(i+1,i+1)=21; 
  end 
  kc(fix(n/2),fix(n/2))=20;
  I = eye(n-2) ; 
  kd = zeros(nnin); %Matriz de nodos internos
  kd(1,1) = -8;
  for i=1:n-3
    kd(i,i+1)=2 ;
    kd(i+1,i)=2 ;
    kd(i+1,i+1)=-8; 
  end 
  k = zeros(nnin*nmin) ; 
  k(1:nnin , 1:nnin) = kb ;
  for j=2:nmin 
    jj = nnin*(j-1)+1;
    k(jj:jj+nnin-1 , jj:jj+nnin-1) = kb;
    k(jj-nnin:jj-1 , jj:jj+nnin-1) = kd;
    k(jj:jj+nnin-1 , jj-nnin:jj-1) = kd;
  end 
  k(nnin+1:nnin+nnin,nnin+1:nnin+nnin)=kc;
  k(1:nnin,nnin*(nnin-1)+1:nnin*nnin)=I;
  k(nnin*(nnin-1)+1:nnin*nnin,1:nnin)=I;


  %-----------------------------------------------------------------------------------------------------------------------
  % Vector de terminos independientes dado por las condiciones de borde y
  % la expresión: (h^4 q_z (x,y))/D, donde h es igual al espaciamiento de
  % los nodos
  f = zeros(nnin*nmin , 1);
  for i = 1:nnin*nmin
    f(i) = (q*hx^4/D); 
%   f(i+nnin*(nmin-1)) = c(i+1 , m) ; 
  end 
  %nuiz = 1:nnin:nnin*nmin ;
  %nude = nnin:nnin:nnin*nmin ; 
  %f(nuiz) = f(nuiz)+c(1,2:m-1)' ;
  %f(nude) = f(nude)+c(n,2:m-1)' ;
%---------------------------------------------------------------------------------------------------------------------------
% Solución del sistema matricial
  u = k\f;

%---------------------------------------------------------------------------------------------------------------------------
% Inclución de los bordes en la solución
  inud = 1;
  for i=2:n-1
    for j=2:m-1
      c(i,j) = u(inud) ;
      inud = inud + 1 ;
    end 
  end 
  
%----------------------------------------------------------------------------------------------------------------------------
% Representación grafica de la placa deflectada.
  x = 0:hx:X ;
  y = 0:hy:Y ;
  colormap("winter");
  surf(x, y, c'),title('Grafica deflexión placa')
  shading("flat")
  a=min(c);
  dmax=min(a)  %Deflexión maxima
  xlabel("Deflexión maxima: "+dmax*1000+ " mm")
  %contour(x, y, ub')
  end 

