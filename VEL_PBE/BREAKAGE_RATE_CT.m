function BREAKAGE = BREAKAGE_RATE_CT(X,Xmin,SIGMA,RHOc,EPSILON,K)

%Coulaloglou & Tavlarides

k1 = 0.336;
k2 = 1.3;

DUMMY = SIGMA/(RHOc*EPSILON^(2/3));

BREAKAGE = k1./ X.^(2/3)*EPSILON^(1.0/3.0) .* exp(- k2*DUMMY./X.^(5/3));

BREAKAGE(X<=Xmin) = 0;

BREAKAGE = K*BREAKAGE;


end%END FUNCTION