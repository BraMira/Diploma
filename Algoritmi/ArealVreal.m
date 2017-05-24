

A = randn(10);
v = randn(10,1)+1i*randn(10,1);
mu = (v'*A*v)/(v'*v);
[n, m] = size(A);
A = A-mu*eye(n);
H=(A+A')/2;
K=(A-A')/(2*1i);
opts.disp = 0; 
%opts.v0 = rand(n,1);
opts.maxit = 1000; 
opts.tol = 10^-4;
napaka = Inf;
korak = 0;

[x, vred_k1] = eigs(K,1,'lr',opts); %najve?ja l. vrednost (realen del) in pripadajo? l. vektor x
[y, vred_k2] = eigs(K,1,'sr',opts);

elipsa(A,x,y);
hold on
[xx, vred_h1] = eigs(H,1,'lr',opts); %najve?ja l. vrednost in pripadajo? l. vektor xx
[yy, vred_h2] = eigs(H,1,'sr',opts);
elipsa(i*A,xx,yy);

hold on
elipsa(A,xx,y);
hold on
elipsa(A,x,yy);
hold off

% A = A + mu*eye(n);
[b nap kor] = izotropniMeurant(A,mu)