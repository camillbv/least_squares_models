function D = LagrangeDerivativeMatrix_GLL(np)

% This method gives the matrix of the lagrange derivatives l'_j(x_i) 
% for i,j=1:np
%
% authors: Kay Hansen-Zahl
%          Tormod Bjo/ntegaard
% 
%          Fall 2002

D = zeros(np,np);
[GLLPoints,GLLWeights]=GaussLobattoLegendre(np);

for i=1:np
  for j=1:np
    if i==j
      D(i,j)=0;
    else
      D(i,j) = LegendreL(np-1,GLLPoints(i))/(LegendreL(np-1,GLLPoints(j))*(GLLPoints(i)-GLLPoints(j)));
    end
  end
end

D(1,1) = -np*(np-1)/4;
D(np,np) = np*(np-1)/4;
    
