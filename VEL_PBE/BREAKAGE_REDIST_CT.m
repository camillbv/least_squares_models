function REDISTb = BREAKAGE_REDIST_CT(x,xp)

Vxp     = pi/6*xp.^3;
P       = pi/2* x^2 * 2.4./Vxp .* exp( -4.5*(2*(pi/6*x^3) - Vxp).^2./Vxp.^2 );
REDISTb = 2*P;

end%END FUNCTION