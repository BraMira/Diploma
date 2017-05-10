function [ b, napaka, korak] = izotropniMeurant(A, mu)
%Èe je možno, ta funkcija izraèuna izotropne vektorje b, ki so rešitev enaèbe mu=b'Ab.
%Input: A...kompleksna matrika
%       mu ... kompleksno število v zalogi vrednosti A
%       
%
%Output: b...izotropni vektor za mu
%        napaka ... norm(b'Ab-mu)
%        korak ... kolikokrat smo raèunali lastne vrednosti in lastne
%        vektorje

if nargin ==1,
    mu=0;
    %tol = 1e-14;
end

[n, m] = size(A);
%preverimo, èe A kvadratna matrika
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
korak = 0;

[x, vred_k1] = eigs(K,1,'lr',opts); %najveèja l. vrednost (realen del) in pripadajoè l. vektor x
[y, vred_k2] = eigs(K,1,'sr',opts); %najmanjša l. vrednost (realen del) in pripadajoè l. vektor y
korak = korak + 1;
pogoj = real(vred_k1)*real(vred_k2);
%èe imamo eno pozitivno in eno negatvno l. vrednost
while (pogoj<0),
    disp('prvi while')
    %poišèemo b1 in b2 s kombiniranjem l. vektorjev
    %kombiniramo s pomoèje xtheta
    [b1, b2] = xtheta(x, y, H, K);
    
    %èe sta nenenièelna jih normiramo
    if sum(abs(b1)<ones(n,1)*1e-10)==0 && sum(abs(b2)<ones(n,1)*1e-10)==0,
        disp('prvi if')
        b1 = b1/norm(b1);
        b2 = b2/norm(b2);
    %uporabimo lemo 3.1.
        if (abs(imag(b1'*A*b1))<1e-10) && (abs(imag(b2'*A*b2))<1e-10),
            disp('drugi if')
            b = lema_31(b1,b2,A);
            if sum(abs(b)<ones(n,1)*1e-10)==0,
                disp('tretji if, konèami pri K')
                napaka = abs(b'*A*b);
                return
            else
                disp('Vektorja b sta enaka 0. Raèunamo s H')
                pogoj = 1;
                continue
            end
        else
            disp('Ne moremo uporabiti leme za vektorje K')
            pogoj = 1;
            continue
        end
    else
        disp('Vektorja sta enaka 0, raèunamo H')
        pogoj = 1;
        continue
    end

end


%ponovimo isti postopek za H in matriko iA
[xx, vred_h1] = eigs(H,1,'lr'); %najveèja l. vrednost in pripadajoè l. vektor xx
[yy, vred_h2] = eigs(H,1,'sr'); %najmanjša l. vrednost in pripadajoè l. vektor yy
korak = korak + 1;
pogoj2 = real(vred_h1)*real(vred_h2);
%èe imamo eno pozitivno in eno negatvno l. vrednost
while (pogoj2<0),
    disp('drugi while')
    %poišèemo b1 in b2 s kombiniranjem l. vrednosti
    %kombiniramo s pomoèje xtheta
    %H = (1i*A +1i*A')/2;
    [b1, b2] = xtheta(xx, yy, H, K);

    %èe sta nenenièelna jih normiramo
    if sum(abs(b1)<ones(n,1)*1e-10)==0 && sum(abs(b2)<ones(n,1)*1e-10)==0,
        disp('prvi if')
        b1 = b1/norm(b1);
        b2 = b2/norm(b2);
    %uporabimo lemo 3.1.
        if (abs(imag(b1'*A*b1))<1e-10) && (abs(imag(b2'*A*b2))<1e-10),
            disp('drugi if')
            b = lema_31(b1,b2,A);
            if sum(abs(b)<ones(n,1)*1e-10)==0,
                disp('tretji if, konèamo pri H')
                napaka = abs(b'*A*b);
                return
            else
                disp('Vektorja b sta enaka 0. Raèunamo s kombinacijo K in H')
                pogoj2 = 1;
                continue
            end
        else
            disp('Ne moremo uporabiti leme za vektorje H')
            pogoj2 = 1;
            continue
        end
    else
        disp('Vektorja sta enaka 0, raèunamo s kombinacijo K in H')
        pogoj2 = 1;
        continue
    end

end

disp('preverja pogoje?')
if (pogoj==1) && (pogoj2==1),
    disp('H ima enako predznaèene lastne vrednosti')
    % uporabimo xthetha za l. vektorja iz H in K
    [b1, b2] = xtheta(x,xx,H,K);
    if sum(abs(b1)<ones(n,1)*1e-10)==0 && sum(abs(b2)<ones(n,1)*1e-10)==0,
        b = [b1, b2];
        napaka = [abs(b1'*A*b1),abs(b2'*A*b2)];
        return
    else
        disp('Funkcija ne najde rešitve.')
    end
end

end
