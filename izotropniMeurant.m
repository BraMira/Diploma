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
K=(A+A')/(2i);
opts.disp = 0; 
opts.maxit = 1000; 
opts.tol = 10^-4;
korak = 0;

[x, vred_k1] = eigs(K,1,'lr',opts); %najveèja l. vrednost in pripadajoè l. vektor x
[y, vred_k2] = eigs(K,1,'sr',opts); %najmanjša l. vrednost in pripadajoè l. vektor y
korak = korak + 1;

%èe so realni deli enaki niè, ne moremo uporabiti algoritma
% if (abs(real(vred_k1)) < tol)==1 || (abs(real(vred_k2))<tol)==1,
%     disp('Algoritma ne moremo uporabiti')
%     return
%èe imamo eno pozitivno in eno negatvno l. vrednost
if real(vred_k1)*real(vred_k2)<0,
    %poišèemo b1 in b2 s kombiniranjem l. vektorjev
    %kombiniramo s pomoèje xtheta
    [b1, b2] = xtheta(x, y, H, K);
    
    %èe sta nenenièelna jih normiramo
%     if (b1~=0)==1 && (b2~=0)==1,
    b1 = b1/norm(b1);
    b2 = b2/norm(b2);
    %uporabimo lemo 3.1.
    if imag(b1'*A*b1)==0 && imag(b2'*A*b2)==0,
        b = lema_31(b1,b2,A);
        napaka = abs(b'*A*b);
    end
%     else
%         return
%     end
else
    disp('K ima enako predznaèene lastne vrednosti')
    %ponovimo isti postopek za H in matriko iA
    [xx, vred_h1] = eigs(H,1,'lr'); %najveèja l. vrednost in pripadajoè l. vektor xx
    [yy, vred_h2] = eigs(H,1,'sr'); %najmanjša l. vrednost in pripadajoè l. vektor yy
    korak = korak + 1;
    %èe imamo eno pozitivno in eno negatvno l. vrednost
    if real(vred_h1)*real(vred_h2)<0,
        %poišèemo b1 in b2 s kombiniranjem l. vrednosti
        %kombiniramo s pomoèje xtheta
        [b1, b2] = xtheta(xx, yy, H, K);

        %èe sta nenenièelna jih normiramo
%         if b1~=0 && b2~=0,
            b1 = b1/norm(b1);
            b2 = b2/norm(b2);
        %uporabimo lemo 3.1.
        if imag(b1'*(1i*A)*b1)==0 && imag(b2'*(1i*A)*b2)==0,
            b = lema_31(b1,b2,1i*A);
            napaka = abs(b'*A*b);
        end
%         else
%             return        
    else
        disp('H ima enako predznaèene lastne vrednosti')
        % uporabimo xthetha za l. vektorja iz H in K
        [b1, b2] = xtheta(x,xx,H,K);
        b = [b1, b2];
        napaka = abs(b'*A*b);
    end
end
end
    
  
    %% Kaj èe predznaka nista ok?
    %%Za koliko l. vektorjev naj to raèunam?


