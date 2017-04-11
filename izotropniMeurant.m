function [ b, napaka, korak] = izotropniMeurant(A, mu)
%�e je mo�no, ta funkcija izra�una izotropne vektorje b, ki so re�itev ena�be mu=b'Ab.
%Input: A...kompleksna matrika
%       mu ... kompleksno �tevilo v zalogi vrednosti A
%       
%
%Output: b...izotropni vektor za mu
%        napaka ... norm(b'Ab-mu)
%        korak ... kolikokrat smo ra�unali lastne vrednosti in lastne
%        vektorje

if nargin ==1,
    mu=0;
    %tol = 1e-14;
end

[n, m] = size(A);
%preverimo, �e A kvadratna matrika
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

[x, vred_k1] = eigs(K,1,'lr',opts); %najve�ja l. vrednost (realen del) in pripadajo� l. vektor x
[y, vred_k2] = eigs(K,1,'sr',opts); %najmanj�a l. vrednost (realen del) in pripadajo� l. vektor y
korak = korak + 1;
pogoj = real(vred_k1)*real(vred_k2);
%�e imamo eno pozitivno in eno negatvno l. vrednost
while (pogoj<0),
    disp('prvi while')
    %poi��emo b1 in b2 s kombiniranjem l. vektorjev
    %kombiniramo s pomo�je xtheta
    [b1, b2] = xtheta(x, y, H, K);
    
    %�e sta neneni�elna jih normiramo
    if (sum(abs(b1)<1e-10)<n) && (sum(abs(b2)<1e-10)<n),
        disp('prvi if')
        b1 = b1/norm(b1);
        b2 = b2/norm(b2);
    %uporabimo lemo 3.1.
        if (abs(imag(b1'*A*b1))<1e-10) && (abs(imag(b2'*A*b2))<1e-10),
            disp('drugi if')
            b = lema_31(b1,b2,A);
            if (sum(abs(b)<1e-10)<n),
                disp('tretji if, kon�ami pri K')
                napaka = abs(b'*A*b);
                return
            else
                disp('Vektorja b sta enaka 0. Ra�unamo s H')
                pogoj = 1;
                continue
            end
        else
            disp('Ne moremo uporabiti leme za vektorje K')
            pogoj = 1;
            continue
        end
    else
        disp('Vektorja sta enaka 0, ra�unamo H')
        pogoj = 1;
        continue
    end

end


%ponovimo isti postopek za H in matriko iA
[xx, vred_h1] = eigs(H,1,'lr'); %najve�ja l. vrednost in pripadajo� l. vektor xx
[yy, vred_h2] = eigs(H,1,'sr'); %najmanj�a l. vrednost in pripadajo� l. vektor yy
korak = korak + 1;
pogoj2 = real(vred_h1)*real(vred_h2);
%�e imamo eno pozitivno in eno negatvno l. vrednost
while (pogoj2<0),
    disp('drugi while')
    %poi��emo b1 in b2 s kombiniranjem l. vrednosti
    %kombiniramo s pomo�je xtheta
    %H = (1i*A +1i*A')/2;
    [b1, b2] = xtheta(xx, yy, H, K);

    %�e sta neneni�elna jih normiramo
    if (sum(abs(b1)<1e-10)<n) && (sum(abs(b2)<1e-10)<n),
        disp('prvi if')
        b1 = b1/norm(b1);
        b2 = b2/norm(b2);
    %uporabimo lemo 3.1.
        if (abs(imag(b1'*A*b1))<1e-10) && (abs(imag(b2'*A*b2))<1e-10),
            disp('drugi if')
            b = lema_31(b1,b2,A);
            if (sum(abs(b)<1e-10)<n),
                disp('tretji if, kon�amo pri H')
                napaka = abs(b'*A*b);
                return
            else
                disp('Vektorja b sta enaka 0. Ra�unamo s kombinacijo K in H')
                pogoj2 = 1;
                continue
            end
        else
            disp('Ne moremo uporabiti leme za vektorje H')
            pogoj2 = 1;
            continue
        end
    else
        disp('Vektorja sta enaka 0, ra�unamo s kombinacijo K in H')
        pogoj2 = 1;
        continue
    end

end

disp('preverja pogoje?')
if (pogoj==1) && (pogoj2==1),
    disp('H ima enako predzna�ene lastne vrednosti')
    % uporabimo xthetha za l. vektorja iz H in K
    [b1, b2] = xtheta(x,xx,H,K);
    if (sum(abs(b1)<1e-10)<n) && (sum(abs(b2)<1e-10)<n)
        b = [b1, b2];
        napaka = [abs(b1'*A*b1),abs(b2'*A*b2)];
        return
    else
        disp('Funkcija ne najde re�itve.')
    end
end

end