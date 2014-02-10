% camilla.berge.vik@ntnu.no, 23.09.2013

%% purpose: model equations for gas plug flow model

function dydz = model_equations(z,y)

%% declare global variables
global eps_G rho_cat v_GS nu a b Mw ptot rho_G factor factor_a factor_b

%% unpack function arguments
w = y;   % weight fractions

%% convert to mole fractions
x = w./Mw/(sum(w./Mw));

%% calculate partial pressures
p = x.*ptot;

%% calculate reaction rate
rxrate = factor*(a*factor_a*p(1)*p(2))/...
    ((1+b*factor_b*p(1))^2);

%% calculate change in mass fractions
dydz = eps_G/v_GS/rho_G*rho_cat.*nu.*rxrate.*Mw;

end