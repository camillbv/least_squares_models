function cEFFICIENCY = COALESCENCE_EFFICIENCY(XI1,XI2,RHOc,EPSILON,SIGMA)

ri  = XI1/2;
rj  = XI2/2;
rij = (1/2)./(1./ri + 1./rj);

CON =  2.302585;
cEFFICIENCY = exp( -CON*RHOc.^0.5.*rij.^(5.0/6.0).*EPSILON^(1.0/3.0)./SIGMA^0.5 );  

end%END FUNCTION