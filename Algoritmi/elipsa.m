function [u,v] = elipsa(x,y,A)
t = linspace(0,2*pi,200);
ko=cos(t);
si=sin(t);
for k=1:200
    z = ko(k).*x+si(k).*y;
    u(k) = real(z'*A*z);
    v(k) = imag(z'*A*z);
end
plot(u,v)
end