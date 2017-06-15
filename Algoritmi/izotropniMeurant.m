function [b, napaka, korak] = izotropniMeurant(A, mu)
%ce je mozno, ta funkcija izracuna izotropne vektorje b, ki so resitev enacbe mu=b'Ab.
%Vhod: A...realna ali kompleksna matrika
%      mu ... kompleksno stevilo v zalogi vrednosti A
%       
%
%Izhod:  b...izotropni vektor za mu
%        napaka ... norm(b'Ab-mu)
%        korak ... kolikokrat smo racunali lastne vrednosti in lastne
%        vektorje
warning off
if nargin ==1,
    mu=0;
    %tol = 1e-14;
end

[n, m] = size(A);
%preverimo, ce A kvadratna matrika
if n~=m,
    disp('Matrika A ni kvadratna!')
    return
end

%problem preobrnemo v b*(A-muI)b=0
A = A-mu*eye(n);

%A realna, algoritem za realne mtr.
if isreal(A)==1,
    [b, napaka, korak] = izotropniMeurantR(A);
    return
end

H=(A+A')/2;
K=(A-A')/(2*1i);
opts.disp = 0; 
opts.maxit = 1000; 
opts.tol = 10^-4;

napaka = Inf;
korak = 0;
b = 0;

[x, vred_k1] = eigs(K,3,'lr',opts); 

[y, vred_k2] = eigs(K,3,'sr',opts); 

korak = korak + 1;

%[X, Y] = vektor(x,y,vred_k1, vred_k2, A);
X = [x,y];
LV = [diag(vred_k1);diag(vred_k2)];
for k=1: size(X,2)-1,
    for j = (k+1):size(X,2),
        if (LV(k,:)*LV(j,:) < 0),
            fi = (log(X(:,k)'*A*X(:,j)*inv(X(:,j)'*A*X(:,k))))/(2i); %#ok<*MINV>
            [b1, b2] = xtheta(X(:,k)*exp(fi*1i), X(:,j), A);


            if sum(abs(b1)<ones(n,1)*1e-10)==0 && sum(abs(b2)<ones(n,1)*1e-10)==0,
                b1 = b1/norm(b1);
                b2 = b2/norm(b2);

                if (abs(imag(b1'*A*b1))<1e-10) && (abs(imag(b2'*A*b2))<1e-10),
                    b = lema_31(b1,b2,A);
                    if sum(abs(b)<ones(n,1)*1e-10)==0,
                        napaka = abs(b'*A*b);
                        return
                    end
                end
            end
        end
    end

end


%ponovimo isti postopek za H in matriko iA
[xx, vred_h1] = eigs(H,3,'lr'); 

[yy, vred_h2] = eigs(H,3,'sr'); 

korak = korak + 1;

XX = [xx,yy];
LV1 = [diag(vred_h1);diag(vred_h2)];

for k=1: size(XX,2)-1,
    for j = (k+1): size(XX,2),
        if (LV1(k,:)*LV1(j,:)< 0),

            fi = (log(XX(:,k)'*A*XX(:,j)*inv(XX(:,j)'*A*XX(:,k))))/(2i); %#ok<*MINV>
            [b1, b2] = xtheta(XX(:,k)*exp(fi*1i), XX(:,j), 1i*A);

            if sum(abs(b1)<ones(n,1)*1e-10)==0 && sum(abs(b2)<ones(n,1)*1e-10)==0,
                b1 = b1/norm(b1);
                b2 = b2/norm(b2);

                if (abs(imag(b1'*(1i*A)*b1))<1e-10) && (abs(imag(b2'*(1i*A)*b2))<1e-10),
                    b = lema_31(b1,b2,(1i*A));
                    if sum(abs(b)<ones(n,1)*1e-10)==0,
                        napaka = abs(b'*(1i*A)*b);
                        return
                    end
                end
            end
        end
    end
end

korak = korak +1;

for k=1:size(X,2),
    for j = 1: size(XX,2),
            fi = (log(X(:,k)'*A*XX(:,j)*inv(XX(:,j)'*A*X(:,k))))/(2i); %#ok<*MINV>
            [b1, b2] = xtheta(X(:,k)*exp(fi*1i), XX(:,j), A);

            if sum(abs(b1)<ones(n,1)*1e-10)==0 && sum(abs(b2)<ones(n,1)*1e-10)==0,
                b1 = b1/norm(b1);
                b2 = b2/norm(b2);

                if (abs(imag(b1'*A*b1))<1e-10) && (abs(imag(b2'*A*b2))<1e-10),
                    b = lema_31(b1,b2,A);
                    if sum(abs(b)<ones(n,1)*1e-10)==0,
                        napaka = abs(b'*A*b);
                        return
                    end
                end
            end
    end
end        
disp('Algoritem ne najde resitve, uporabi algoritem CPU')

%izboljsava
z1 = [real(x(:,1)'*A*x(:,1)), imag(x(:,1)'*A*x(:,1))];
z2 = [real(y(:,1)'*A*y(:,1)), imag(y(:,1)'*A*y(:,1))];
%y = a*x+c
a = (z2(2)-z1(2))/(z2(1)-z1(1));
c = z1(2) - z1(1)*a;
nicla = - c/a;
if nicla <0,
    z3 = [real(xx(:,1)'*A*xx(:,1)), imag(xx(:,1)'*A*xx(:,1))];
    a1 = (z3(2)-z1(2))/(z3(1)-z1(1));
    c1 = z3(2) - z3(1)*a1;
    nicla1 = - c1/a1;
    a2 = (z3(2)-z2(2))/(z3(1)-z2(1));
    c2 = z3(2) - z3(1)*a2;
    nicla2 = - c2/a2;
    if nicla1 >0 && nicla1 <= z3(1),
%         bb = @(t,th) exp(-1i*th)*nicla + t*nicla1; %b(t,theta)
%         al = @(th) exp(1i*th)*nicla'*A*nicla1 + exp(-1i*th)*nicla1'*A*nicla; %alfa(theta)
% 
%         th = angle(nicla1'*A*nicla - nicla.'*conj(A)*conj(nicla1));
%         t1 = (-al(th) + sqrt(al(th)^2 -4*alfa1*alfa2))/(2*alfa2);
% 
%         b = bb(t1,th)/norm(bb(t1,th));
        b = lema_31(
        napaka = abs(b'*A*b);
        return
    elseif nicla2>0 && nicla2<=z3(1),
        b = lema_31(nicla,nicla2,A);
        napaka = abs(b'*A*b);
        %return
    else
        disp('Uporabi algoritem Cardna')
    end
elseif nicla >0
    z3 = [real(yy(:,1)'*A*yy(:,1)), imag(yy(:,1)'*A*yy(:,1))];
    a1 = (z3(2)-z1(2))/(z3(1)-z1(1));
    c1 = z3(2) - z3(1)*a1;
    nicla1 = - c1/a1;
    a2 = (z3(2)-z2(2))/(z3(1)-z2(1));
    c2 = z3(2) - z3(1)*a2;
    nicla2 = - c2/a2;
    if nicla1 <0 && nicla1 >= z3(1),
        b = lema_31(nicla/norm(nicla),nicla1/norm(nicla1),A);
        napaka = abs(b'*A*b);
        %return
    elseif nicla2<0 && nicla2 >= z3(1),
        b = lema_31(nicla/norm(nicla),nicla2/norm(nicla2),A);
        napaka = abs(b'*A*b);
        %return
    else
        disp('Uporabi algoritem Cardna')
    end
end
    
        



end
