A  = randn(5) + i*randn(5);
H=(A+A')/2;
K=(A-A')/(2*1i);
opts.disp = 0; 
opts.maxit = 1000; 
opts.tol = 10^-4;
korak = 0;

[x, vred_k1] = eigs(K,1,'lr',opts); %najveèja l. vrednost (realen del) in pripadajoè l. vektor x
[y, vred_k2] = eigs(K,1,'sr',opts); %najmanjša l. vrednost (realen del) in pripadajoè l. vektor y


t = linspace(0,2*pi,200);
ko=cos(t);
si=sin(t);
for k=1:200
    z = ko(k).*x+si(k).*y;
    u(k) = real(z'*A*z);
    v(k) = imag(z'*A*z);
end

plot(u,v)