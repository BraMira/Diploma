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


for k=1: size(X,2),
    for j = 1:size(X,2),
        if k~=j,
        
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
                    else
                        disp('Vektorja b sta enaka 0. Racunamo s H')
                        continue
                    end
                else
                    disp('Ne moremo uporabiti leme za vektorje K')
                    continue
                end
            else
                disp('Vektorja sta enaka 0, racunamo H')
                continue
            end
        end
    end

end


%ponovimo isti postopek za H in matriko iA
[xx, vred_h1] = eigs(H,3,'lr'); 

[yy, vred_h2] = eigs(H,3,'sr'); 

korak = korak + 1;

%[XX, YY] = vektor(xx,yy,vred_h1, vred_h2, 1i*A);
XX = [xx,yy];

for k=1: size(XX,2),

    for j = 1: size(XX,2),
        if k~=j,
            

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
                    else
                        disp('Vektorja b sta enaka 0. Racunamo s kombinacijo K in H')
                        continue
                    end
                else
                    disp('Ne moremo uporabiti leme za vektorje H')
                    continue
                end
            else
                disp('Vektorja sta enaka 0, racunamo s kombinacijo K in H')
                continue
            end
        end
    end
end


XXX = [x,y];
YYY = [xx,yy];


korak = korak +1;

for k=1:size(XXX,2),
    for j = 1: size(YYY,2),
        fi = (log(XXX(:,k)'*A*YYY(:,j)*inv(YYY(:,j)'*A*XXX(:,k))))/(2i); %#ok<*MINV>
        [b1, b2] = xtheta(XXX(:,k)*exp(fi*1i), YYY(:,j), A);

        if sum(abs(b1)<ones(n,1)*1e-10)==0 && sum(abs(b2)<ones(n,1)*1e-10)==0,
            b1 = b1/norm(b1);
            b2 = b2/norm(b2);

            if (abs(imag(b1'*A*b1))<1e-10) && (abs(imag(b2'*A*b2))<1e-10),
                b = lema_31(b1,b2,A);
                if sum(abs(b)<ones(n,1)*1e-10)==0,
                    napaka = abs(b'*A*b);
                    return
                else
                    disp('Vektorja b sta enaka 0')
                    continue
                end
            else
                disp('Ne moremo uporabiti leme')
                continue
            end
        else
            disp('Vektorja sta enaka 0')
            continue
        end
    end
end        
disp('Algoritem ne najde rešitve, uporabi algoritem CPU')

end
