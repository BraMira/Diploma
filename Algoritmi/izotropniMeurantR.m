function [b, napaka, korak] = izotropniMeurantR(A)
%Resi problem iztotropnih vektorjev, èe je matrika A realna in mi realen
%oz. 0

%racunamo b'Hb=0
H=(A+A')/2;

opts.disp = 0; 
opts.maxit = 1000; 
opts.tol = 10^-4;

[x, vred_h1] = eigs(H,1,'la',opts); 

[y, vred_h2] = eigs(H,1,'sa',opts); 

if (vred_h1 * vred_h2) < 0,
    b1 = sqrt(vred_h1 / (vred_h1 + abs(vred_h2)))*y + sqrt( abs(vred_h2) / (vred_h1 + abs(vred_h2)))*x;
    b2 = -sqrt(vred_h1 / (vred_h1 + abs(vred_h2)))*y + sqrt( abs(vred_h2) / (vred_h1 + abs(vred_h2)))*x;
    b = [b1 ,b2];
    napaka = [abs(b1'*H*b1),abs(b2'*H*b2)];
    korak = 1;
else
    b = 0;
    napaka = Inf;
    korak = 0;
    disp('H je definitna')
end

