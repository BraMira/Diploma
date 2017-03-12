function [ b, napaka, korak] = izotropniMeurant(A, mu )
%Èe je možno, ta funkcija izraèuna izotropne vektorje b, ki so rešitev enaèbe mu=b'Ab.
%Input: A...kompleksna matrika
%       mu ... kompleksno število v zalogi vrednosti A
%
%Output: b...izotropni vektor za mu
%        napaka ... norm(b'Ab-mu)
%        korak ... kolikokrat smo raèunali lastne vrednosti in lasnte
%        vektorje

if nargin ==1,
    mu=0;
end

[n m] = size(A);
%preverimo, èe A kvadratna matrika
if n~=m,
    disp('Matrika A ni kvadratna!')
    return
end

%problem preobrnemo v b*(A-muI)b=0
A = A-mu*eye(n);

H=(A+A')/2;
K=(A+A')/(2i);

korak = 0;

[x, vred_k1] = eigs(K,1,'lr'); %najveèja l. vrednost in pripadajoè l. vektor x
[y, vred_k2] = eigs(K,1,'sr'); %najmanjša l. vrednost in pripadajoè l. vektor y
korak = korak + 1;

if sign(real(vred_k1))~=sign(real(vred_k2)),
    %disp('sem na prvem koraku')
%%  Preverim èe sta lastni vrednosti razl. predznaka?

    %Kombinacija vektorjev K
    alfa = imag(x'*H*x + i*x'*K*x);
    beta = imag(y'*H*y + i*y'*K*y);
    gama = imag(x'*H*y + y'*H*x + i*(x'*K*y + y'*K*x));

    %rešimo enaèbo beta*t^2 + gama*t + alfa = 0
    t1 = (-gama + sqrt(gama^2 -4*beta*alfa))/(2*beta);
    t2 = (-gama - sqrt(gama^2 -4*beta*alfa))/(2*beta);

    %%  preverim èe real(t1)>0?

    %t = tan(theta) (predpostavili smo, da cos(theta)~=0
    theta1 = atan(t1);
    theta2 = atan(t2);
    %%potrebno preveriti èe je theta med 0 in pi?

    b1 = cos(theta1)*x + sin(theta1)*y;
    b2 = cos(theta2)*x + sin(theta2)*y;
    
    %% ali sta b1, b2 enotska? Moram jaz normirat ali morata tak že biti
    b1 = b1/norm(b1);
    b2 = b2/norm(b2);
    alfa1 = real(b1'*A*b1);
    alfa2 = real(b2'*A*b2);
    
    if sign(alfa1)==1 && sign(alfa2)==-1,
        alfa1 = alfa2;
        alfa2 = real(b1'*A*b1);
    end
    %% Kaj èe predznaka nista ok?
    %% Kaj, èe imaginarna dela b'Ab nista 0?
    if imag(b1'*A*b1) ~=0,
        disp('ni 0')
    end
    if imag(b2'*A*b2)~=0,
        disp('b2 ni 0')
    end
    %% Uporabimo lemo 3.1.
    
    bb = @(t,th) exp(-i*th)*b1 + t*b2;
    al = @(th) exp(i*th)*b1'*A*b2 + exp(-i*th)*b2'*A*b1;
    
    th = angle(b2'*A*b1 + b1.'*conj(A)*conj(b2));
    tt1 = (-al(th) + sqrt(al(th)^2 -4*alfa1*alfa2))/(2*alfa2);
    
    b = bb(tt1,th)/norm(bb(tt1,th));

    napaka = abs(b'*A*b);

else
    disp('K nima razlièno predznaèenih lastnih vrednosti')
    disp('Preverimo lastne vrednosti H in delamo z matriko iA')
    [x, vred_k1] = eigs(H,1,'lr'); %najveèja l. vrednost in pripadajoè l. vektor x
    [y, vred_k2] = eigs(H,1,'sr'); %najmanjša l. vrednost in pripadajoè l. vektor y
    A = i*A;
    korak = korak + 1;

    if sign(real(vred_k1))~=sign(real(vred_k2)),
    %  Preverim èe sta lastni vrednosti razl. predznaka?
        %Kombinacija vektorjev 
        alfa = imag(x'*H*x + i*x'*K*x);
        beta = imag(y'*H*y + i*y'*K*y);
        gama = imag(x'*H*y + y'*H*x + i*(x'*K*y + y'*K*x));

        %rešimo enaèbo beta*t^2 + gama*t + alfa = 0
        t1 = (-gama + sqrt(gama^2 -4*beta*alfa))/(2*beta);
        t2 = (-gama - sqrt(gama^2 -4*beta*alfa))/(2*beta);

        %%  preverim èe real(t1)>0?

        %t = tan(theta) (predpostavili smo, da cos(theta)~=0
        theta1 = atan(t1);
        theta2 = atan(t2);
        %%potrebno preveriti èe je theta med 0 in pi?

        b1 = cos(theta1)*x + sin(theta1)*y;
        b2 = cos(theta2)*x + sin(theta2)*y;

        %% ali sta b1, b2 enotska? Moram jaz normirat ali morata tak že biti
        b1 = b1/norm(b1);
        b2 = b2/norm(b2);
        alfa1 = real(b1'*A*b1);
        alfa2 = real(b2'*A*b2);

        if sign(alfa1)==1 && sign(alfa2)==-1,
            alfa1 = alfa2;
            alfa2 = real(b1'*A*b1);
        end
        %% Kaj èe predznaka nista ok?
        %% Kaj, èe imaginarna dela b'Ab nista 0?
        if imag(b1'*A*b1) ~=0,
            disp('ni 0')
        end
        if imag(b2'*A*b2)~=0,
            disp('b2 ni 0')
        end
        %% Uporabimo lemo 3.1.

        bb = @(t,th) exp(-i*th)*b1 + t*b2;
        al = @(th) exp(i*th)*b1'*A*b2 + exp(-i*th)*b2'*A*b1;

        th = angle(b2'*A*b1 + b1.'*conj(A)*conj(b2));
        tt1 = (-al(th) + sqrt(al(th)^2 -4*alfa1*alfa2))/(2*alfa2);

        b = bb(tt1,th)/norm(bb(tt1,th));

        napaka = abs(b'*A*b);
    end
    %%else èe niti ta konstrukcija ne deluje

end

end

