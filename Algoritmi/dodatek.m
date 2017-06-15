
Vbox = [x(:,1), xx(:,1), y(:,1),yy(:,1)];
rvbox = [x(:,1)'*A*x(:,1), xx(:,1)'*A*xx(:,1), y(:,1)'*A*y(:,1), yy(:,1)'*A*yy(:,1)];
tol = 1e-14;

b = izboljsava(A,Vbox,rvbox,tol);
napaka = abs(b'*A*b);
korak = 4;