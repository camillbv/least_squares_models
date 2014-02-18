function cFREQUENCY  = COALESCENCE_FREQUENCY(XI1,XI2,RHOc,EPSILON,SIGMA,K)

cRATE       = COALESCENCE_RATE(XI1,XI2,EPSILON);
cEFFICIENCY = COALESCENCE_EFFICIENCY(XI1,XI2,RHOc,EPSILON,SIGMA);

cFREQUENCY  = cRATE.*cEFFICIENCY;

cFREQUENCY  = K*cFREQUENCY;

end%END FUNCTION