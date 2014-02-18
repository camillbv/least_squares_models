function h = BREAKAGE_REDIST_DIEMER_OLSON(xi,pxi)

    xi  = xi./pxi;

    nu = 2;
    q  = 2;

    P = 3*gamma(q*nu)/gamma(q)/gamma(q*(nu-1))*(xi.^3).^(q-1).*(1-xi.^3).^(q*(nu-1)-1).*xi.^2;
    P = P./pxi;
    
    h = P*nu;
end