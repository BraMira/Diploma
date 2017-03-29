function [b1,b2] = xtheta(x,y,H,K)
%
alfa = imag(x'*H*x + 1i*x'*K*x);
beta = imag(y'*H*y + 1i*y'*K*y);
gama = imag(x'*H*y + y'*H*x + 1i*(x'*K*y + y'*K*x));

%re�imo ena�bo beta*t^2 + gama*t + alfa = 0
t1 = (-gama + sqrt(gama^2 -4*beta*alfa))/(2*beta);
t2 = (-gama - sqrt(gama^2 -4*beta*alfa))/(2*beta);

%preverimo, �e sta re�itvi realni
if (abs(imag(t1))<1e-10),
    theta1 = atan(real(t1));%t = tan(theta) (predpostavili smo, da cos(theta)~=0
    b1 = cos(theta1)*x + sin(theta1)*y;
else
    b1 = 0;
end

if (abs(imag(t2))<1e-10),
    theta2 = atan(real(t2));
    b2 = cos(theta2)*x + sin(theta2)*y;
else
    b2 = 0;
end
%%potrebno preveriti �e je theta med 0 in pi?
end