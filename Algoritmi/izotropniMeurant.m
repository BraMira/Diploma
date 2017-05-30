function [ b, napaka, korak] = izotropniMeurant(A, mu)
%ce je mozno, ta funkcija izracuna izotropne vektorje b, ki so resitev enacbe mu=b'Ab.
%Input: A...kompleksna matrika
%       mu ... kompleksno stevilo v zalogi vrednosti A
%       
%
%Output: b...izotropni vektor za mu
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

H=(A+A')/2;
K=(A-A')/(2*1i);
opts.disp = 0; 
opts.maxit = 1000; 
opts.tol = 10^-4;
napaka = Inf;
korak = 0;

[x, vred_k1] = eigs(K,1,'lr',opts); 
%najvecja l. vrednost (realen del) in pripadajoc l. vektor x
[y, vred_k2] = eigs(K,1,'sr',opts); 
%najmanjsa l. vrednost (realen del) in pripadajoc l. vektor y
korak = korak + 1;
pogoj = real(vred_k1)*real(vred_k2);
%ce imamo eno pozitivno in eno negatvno l. vrednost
while (pogoj<0),
    %poiscemo b1 in b2 s kombiniranjem l. vektorjev
    %kombiniramo s pomocje xtheta
    [b1, b2] = xtheta(x, y, A);
    
    %ce sta nenenicelna jih normiramo
    if sum(abs(b1)<ones(n,1)*1e-10)==0 && sum(abs(b2)<ones(n,1)*1e-10)==0,
        b1 = b1/norm(b1);
        b2 = b2/norm(b2);
    %uporabimo lemo 3.1.
        if (abs(imag(b1'*A*b1))<1e-10) && (abs(imag(b2'*A*b2))<1e-10),
            b = lema_31(b1,b2,A);
            if sum(abs(b)<ones(n,1)*1e-10)==0,
                napaka = abs(b'*A*b);
                return
            else
                disp('Vektorja b sta enaka 0. Racunamo s H')
                pogoj = 1;
                continue
            end
        else
            disp('Ne moremo uporabiti leme za vektorje K')
            pogoj = 1;
            continue
        end
    else
        disp('Vektorja sta enaka 0, racunamo H')
        pogoj = 1;
        continue
    end

end


%ponovimo isti postopek za H in matriko iA
[xx, vred_h1] = eigs(H,1,'lr'); %najvecja l. vrednost in pripadajoc l. vektor xx
[yy, vred_h2] = eigs(H,1,'sr'); %najmanjsa l. vrednost in pripadajoc l. vektor yy
korak = korak + 1;
pogoj2 = real(vred_h1)*real(vred_h2);
%ce imamo eno pozitivno in eno negatvno l. vrednost
while (pogoj2<0),
    %poiscemo b1 in b2 s kombiniranjem l. vrednosti
    %kombiniramo s pomocje xtheta
    %H = (1i*A +1i*A')/2;
    [b1, b2] = xtheta(xx, yy, (1i*A));

    %ce sta nenenicelna jih normiramo
    if sum(abs(b1)<ones(n,1)*1e-10)==0 && sum(abs(b2)<ones(n,1)*1e-10)==0,
        b1 = b1/norm(b1);
        b2 = b2/norm(b2);
    %uporabimo lemo 3.1.
        if (abs(imag(b1'*(1i*A)*b1))<1e-10) && (abs(imag(b2'*(1i*A)*b2))<1e-10),
            b = lema_31(b1,b2,(1i*A));
            if sum(abs(b)<ones(n,1)*1e-10)==0,
                napaka = abs(b'*(1i*A)*b);
                return
            else
                disp('Vektorja b sta enaka 0. Racunamo s kombinacijo K in H')
                pogoj2 = 1;
                continue
            end
        else
            disp('Ne moremo uporabiti leme za vektorje H')
            pogoj2 = 1;
            continue
        end
    else
        disp('Vektorja sta enaka 0, racunamo s kombinacijo K in H')
        pogoj2 = 1;
        continue
    end

end
disp('preverimo pogoje')
if (pogoj==1) && (pogoj2==1),
    korak = korak +1;
    % uporabimo xthetha za l. vektorja iz H in K
    [b1, b2] = xtheta(x,yy,A);
    disp('ali vektorja enaka 0?')
    if sum(abs(b1)<ones(n,1)*1e-10)==0 && sum(abs(b2)<ones(n,1)*1e-10)==0,
        
        b1 = b1/norm(b1);
        b2 = b2/norm(b2);
        disp('ali produkt realen?')
        if (abs(imag(b1'*A*b1))<1e-10) && (abs(imag(b2'*A*b2))<1e-10),
            
            b = lema_31(b1,b2,A);
            disp('ali nismo ni?elen vektor dobili?')
            if sum(abs(b)<ones(n,1)*1e-10)==0,
                napaka = abs(b'*A*b);
                return
            else
                disp('Funkcija ne najde resitve.')
            end
        else
            disp('Funkcija ne najde resitve.')
        end
        
%     if sum(abs(b1)<ones(n,1)*1e-10)==0 && sum(abs(b2)<ones(n,1)*1e-10)==0,
%         b = [b1, b2];
%         napaka = [abs(b1'*A*b1),abs(b2'*A*b2)];
        %return
    else
        disp('Funkcija ne najde resitve.')
    end
end

end
