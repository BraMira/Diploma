plot(0,0, '*')
plot(x(:,1)'*A*x(:,1),'*r')
plot(x(:,2)'*A*x(:,2),'*r')
plot(x(:,3)'*A*x(:,3),'*r')
plot(y(:,1)'*A*y(:,1),'*g')
plot(y(:,2)'*A*y(:,2),'*g')
plot(y(:,3)'*A*y(:,3),'*g')

%elipsa...
[b1, b2] = xtheta(x(:,1), y(:,2), A); %ki sta najbolj levo
 plot(b1'*A*b1,0,'+r')
 plot(b2'*A*b2,0,'+r')
 b1 = b1/norm(b1);
 b2 = b2/norm(b2);
 b = lema_31(b1,b2,A)
 
 plot(b'*A*b, '+b')