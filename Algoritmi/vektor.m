function [x,y] = vektor(x_maxlv,y_minlv,lv1,lv2,A)
%Funkcija vektor vrne tiste vektorje x,y, ki so najbolj primerni za
%konstrukcijo elipse
%       x_maxlv.....najvecji lastni vektorji
%       y_minlv.....najmanjsi lastni vektorji
%       lv1.....lastne vrednosti pripadajoce najvecjim lastnim vektorjem
%       lv2.....lastne vrednosti pripadajoce najmanjsim lastnim vektorjem

%lastne vrednosti morajo biti nasprotno predznaèene
x1 = x_maxlv(:,diag(lv1>0));
y1= y_minlv(:,diag(lv2<0));

z1 = diag(x1'*A*x1);
z2 = diag(y1'*A*y1);

%tocke morajo biti zgoraj in spodaj realne osi
x1 = x1(:, imag(z1)>0);
z1 = z1(imag(z1)>0);

y1 = y1(:, imag(z2)<0);
z2 = z2(imag(z2)<0);

%iscemo tocke najblizje nicli, glede na realno os
min_x = min(abs(real(z1)));
min_y = min(abs(real(z2)));

%vrne vektorje, ki so najblizje 0 glede na realno os in tudi vektorje, ki
%so od njih oddaljeni za manj kot 0.2 glede na realno os
x = x1(:, abs( abs(real(z1)) - min_x) <0.2);
y = y1(:, abs( abs(real(z2)) - min_y) <0.2);
end

    
        





