function Y=nearritz(H,z,tol)
%Determines the vector Y such that z=Y'*HY
%H must be a 2x2 matrix
%tol tolerance for determining if H is normal
%Uses Horn and Johnson analysis of field of values of 2x2 matrix from
%Topics in Matrix Analysis.
%Things decouple nicely for 2x2 matrices of the form
%[ 0 a ; b 0 ]
[SH]=eig(H);
cH=trace(H)/2;
tH=SH-cH;
rangle=angle(tH(1,1));
%save thedata4
tH=(H-cH*eye(2))*exp(-1i*rangle);
[U S]=schur(tH);
CHm= S;
Dm=diag([1,exp(-1i*angle(CHm(1,2)))]);
CH=(Dm'*CHm*Dm);
CH=real(CH);

a=CH(2,2);
b=CH(1,2)/2;
d=sqrt(a^2+b^2);
if abs(a)<tol &&abs(b)<tol
    %Matrix has only one normal eigenvalue
    Y=randn(2,1);
    Y=Y/norm(Y);
    return
end