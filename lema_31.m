function [b] = lema_31(b1,b2,A)
%Input: enotska vektorja b1 in b2, za katera velja imag(bi'*A*bi)=0 za
%i=1,2, za dano matriko A
%Output: izotropni vektor b

if real(b1'*A*b1)<0,
    if real(b2'*A*b2)>0,
        alfa1 = real(b1'*A*b1);
        alfa2 = real(b2'*A*b2);
    else
        disp('b1 in b2 imata enako (negativno) predznaèena realna dela')
        b=0;
        return
    end
elseif real(b1'*A*b1)>0,
    if real(b2'*A*b2)<0,
        alfa2 = real(b1'*A*b1);
        alfa1 = real(b2'*A*b2);
    else
        disp('b1 in b2 imata enako (pozitivno) predznaèena realna dela, ne moremo ju uporabiti!')
        b=0;
        return
    end
else
    disp('ne moremo uporabiti leme')
    b = 0;
    return
end
bb = @(t,th) exp(-1i*th)*b1 + t*b2; %b(t,theta)
al = @(th) exp(1i*th)*b1'*A*b2 + exp(-1i*th)*b2'*A*b1; %alfa(theta)

th = angle(b2'*A*b1 + b1.'*conj(A)*conj(b2));
t1 = (-al(th) + sqrt(al(th)^2 -4*alfa1*alfa2))/(2*alfa2);

b = bb(t1,th)/norm(bb(t1,th));

end