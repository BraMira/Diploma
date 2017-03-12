function [vf eeval]=inversefov(A,z,use_eigs,tol,iters)
%If possible, this function determines a solution to the
%inverse field of values problem: Given a z in the field
%of values of A determine a vector vf such that
%z = vf'*A*vf.
%  Inputs:   A a square matrix
%            z a complex number in the field of values of A
%     use_eigs 1 indicates whether eigs will be used rather than eig
%          tol tolerance used for deciding how  small |vf'*A*vf-z|
%              must be for accepting vf as a generating vector for z
%        iters Maximum number of iterations
%
%  Outputs: vf generating vector for z
%        eeval Number of eigenvalue computations performed
%
%
% R. Carden February 6 2011
%


warning off

if nargin<2
    disp('Must specify a matrix A and a point in the numerical range of A')
    vf=[];
    eeval=0;
    
    return
elseif nargin <3
    tol=1e-14;
    use_eigs=0;    
elseif nargin <4
    tol=1e-14;
elseif nargin <5
    iters=30;
end

if use_eigs==1
    boundarypt = @eigs_boundarypt;
else
    boundarypt = @eig_boundarypt;
end

opts.disp=0;
opts.tol=tol/1000;
err1=1000;
err0=1;
[ n m]=size(A);
jj=1;
eeval=0; %Number of eigenvalue evaluations performed
Vbox=zeros(n,4);
rvbox=zeros(1,4);

%Determine four points on the boundary of numerical range
for i=1:2
    th=((i-1)*pi/2)-pi;
    AA=exp(-1i*th)*(A);
    HA=(AA+AA')/2;
    if use_eigs==0
        [va sa]=eig(HA);eeval=eeval+1;
        sa=diag(sa);
    else
        [v1 s1]=eigs(HA,1,'LR',opts);eeval=eeval+1;
        [v2 s2]=eigs(HA,1,'SR',opts);eeval=eeval+1;
        va=[v1 v2];
        sa=[s1 s2];
    end
    [val I ]= max(sa);
    Vbox(:,1+(i-1))=va(:,I);
    rvbox(1+(i-1))=va(:,I)'*A*va(:,I);
    [val I ]= min(sa);
    Vbox(:,3+(i-1))=va(:,I);
    rvbox(3+(i-1))=va(:,I)'*A*va(:,I);
end


%pick a point outside the numerical range for determining if $z$ is in
%polygonal approximation
q=real(rvbox(3))+ max(.1,.1*real(rvbox(3)-rvbox(1)))+1i*(imag(rvbox(4))+ max(.1,.1*imag((rvbox(4)-rvbox(2)))));

%Form outer approximation to numerical range
x1=min(real(rvbox));y1=min(imag(rvbox));
x2=max(real(rvbox));y2=max(imag(rvbox));
obox=[x1+1i*y1 x2+1i*y1 x2+1i*y2 x1+1i*y2];

ttol=max(x2-x1,y2-y1)*tol*2;
if  ((x2-x1)<=1e-14*(x2+x1)/2) && ((y2-y1)<=1e-14*(y2+y1)/2)
    ttol=tol;
end
tol=ttol;
%Check if matrix is skew-symmetric, symmetric or if numerical range
%is a point.
if ((x2-x1)<tol ||(y2-y1)<tol)
    %If the outer approximation is just a line then
    %Then we need only two points to determine the numerical range
    %First eleiminate unnecessary points
    if ~((x2-x1)<1e-14*(x2+x1)/2 &&(y2-y1)<(1e-14*(y2+y1)/2))
        D=-ones(4);
        for i =1:3
            for j=i+1:4
                D(i,j)=abs(rvbox(i)-rvbox(j));
            end
        end
        [val J ] = max(D);
        [val2 I] = max(val);
        Vc=orth(Vbox(:,[I J(I)]));
        
        H=Vc'*A*Vc;
        
        Y=nearritz(H,z,tol);
        vf=Vc*Y;
        rvf=vf'*A*vf;
    else   %Numerical Range is a point
        rvbox=rvbox(1);
        Vbox=Vbox(:,1);
        rvf=rvbox;
        vf=Vbox;
    end
    
    %If numerical range is a line determine closest point to z
    if abs(rvf-z)>tol
        rvf=[];
        vf=[];
        eeval=-eeval;
        %   disp('Clearly not in the numerical range')
    end
else
    %If we are here then the area of the outer approximation is not
    %negligable.  Check that z is inside the outer approximation.
    [ni fon] =insidequad([x1+1i*y1 x2+1i*y1 x2+1i*y2 x1+1i*y2],z,q,tol);
    if mod(ni,2)==1
        ni=0;
        %Determine if inner approximation is a line
        drvbox=[rvbox rvbox(1)];
        drvbox=abs(diff(drvbox));
        i=find(drvbox<tol);
        rvbox(i)=[];
        Vbox(:,i)=[];
        if length(rvbox)>2
            [ni fon] =insidequad(rvbox,z,q,tol);
        end
        if mod(ni,2)==1
            err1=0;
            err0=0;
            %z is strictly inside the inner approximation
            %Problem reduces to determining a boundary point rvc of the inner
            %approximation such that rvbox(1)-rvc passes through z
            angle1=angle((z-rvbox(1))*conj(rvbox(3)-rvbox(1)));
            
            
            if abs(angle1)<1e-14   %If z lies on rvbox(1)-rvbox(3)
                Vc=orth(Vbox(:,[1 3]));
                rvc=z;
            else
                if angle1>0 &&length(rvbox)~=3
                    rvbox=rvbox([1 3 4]);
                    Vbox=Vbox(:,[1 3 4]);
                end
                %Use law of sines to determine where ray from rvbox to z intersects.
                %ray from rvbox(3)to rvbox(2)
                angle1=angle((rvbox(3)-rvbox(1))*conj(z-rvbox(1)));
                angle2=angle((rvbox(2)-rvbox(3))*conj(rvbox(1)-rvbox(3)));
                angle3=pi-angle1-angle2;
                
                rv3torvc=sin(angle1)/(sin(angle3)/abs(rvbox(3)-rvbox(1)));
                rvc=rv3torvc*exp(1i*angle(rvbox(2)-rvbox(3)))+rvbox(3);
                Vc=orth(Vbox(:,[2 3]));
            end
            %Determine Ritz vector for rvc
            
            H=Vc'*A*Vc;
            Y=nearritz(H,rvc,tol);
            vc=Vc*Y;
            rvc=vc'*A*vc;
            
            %Determine a Ritz vector for z
            if abs(rvc-z)<tol
                vf=vc;
                rvf=rvc;
            else
                
                Vc=orth([vc Vbox(:,1)]);
                H=Vc'*A*Vc;
                Y=nearritz(H,z,tol);
                vf=Vc*Y;
                rvf=vf'*A*vf;
            end
        else
            
            %First determine which side of inner approximation to which z
            % is closest.Must center about trace to prevent undefined angles.
            trA=trace(A)/n;
            theangles=myunwrap(angle(rvbox-trA));
            zangle= myunwrap([angle(rvbox(1)-trA) angle(z-trA)]);
            zangle=zangle(2);
            theangles = theangles-zangle;
            an=sum(theangles<0);
            if an==0
                an=length(theangles);
            end
            bn=mod(an,length(theangles))+1;
            %Denote vertices of new triangle in numerical range
            %with rva rvb rvc
            %Where rva and rvb are determine by vertices of side from inner
            %approximation.  rvc is determined below.
            %rvf will be a ritz value close to z
            
            rva=rvbox(an);
            rvb=rvbox(bn);
            va=Vbox(:,an);
            vb=Vbox(:,bn);
            rvf=rva;
            err1=abs(z-rvf);
            while abs(z-rvf)>tol
                %Determine closest point on rvarvb
                rvd=real((z-rva)*conj(rvb-rva))/abs(rvb-rva);
                rvd=rvd*exp(1i*angle(rvb-rva))+rva;
                rvf=rvd;
                %The angle from x1 to z must always gives us a new point on the
                %boundary of the numerical range unless of
                %course z is on the line between rva and rvb
                %Might want to check if z and rvd are indeed the same point
                %because the angle would be undefined in that case
                if abs(z-rvf)>tol
                    th=angle(z-rvd);
                    AA=exp(-1i*th)*(A);
                    HA=(AA+AA')/2;
                    vc=boundarypt();
                    rvc=vc'*A*vc;
                    eeval=eeval+1;
                    
                    %Now we have a new point and we can immediately
                    %check whether z is in this new triangle or
                    %if we must select a new base for our next triangle.
                    %Outer approximation check
                    if  tol<real(z*exp(-1i*th))-(real(rvc*exp(-1i*th)))
                        rvf=[];
                        vf=[];
                        eeval=-eeval;
                        break;
                    end
                    
                    [ni fon] =insidequad([rva rvb rvc],z,q,tol);
                    if mod(ni,2)==1
                        %Use law of sines to determine x1, where ray from rva to z
                        %intersects opposing side of triangle.
                        %If the triangle is skinny this problem may
                        %be ill-conditioned, in either case we keep
                        %only the best approximation to z.
                        
                        angle1=angle((rvb-rva)*conj(z-rva));
                        angle2=angle((rvc-rvb)*conj(rva-rvb));
                        angle3=pi-angle1-angle2;
                        rvbtorvd=sin(angle1)/(sin(angle3)/abs(rvb-rva));
                        rvd=rvbtorvd*exp(1i*angle(rvc-rvb))+rvb;
                        %Determine Ritz vector for rvd
                        Vc=orth([vb vc]);
                        H=Vc'*A*Vc;
                        Y=nearritz(H,rvd,tol);
                        vd=Vc*Y;
                        rvd=vd'*A*vd;
                        if abs(rvd-z)<tol
                            vf=vd;
                            rvf=rvd;
                            break
                        else
                            %Determine a Ritz vector for this point z
                            Vc=orth([va vd]);
                            H=Vc'*A*Vc;
                            Y=nearritz(H,z,tol);
                            vf=Vc*Y;
                            rvf=vf'*A*vf;
                            if abs(rvf-z)>abs(rvd-z)
                                %In this case we encountered an ill
                                %conditioned problem (skinny triangle) in determining rvd
                                %because rvf should have exactly generated
                                %z.
                                vf=vd;
                                rvf=rvd;
                            end
                            break
                        end
                    end
                    
                else
                    %Determine Ritz vector for point nearest to z
                    %This is the only part of the code that if possible takes
                    %advantage of nondegenerate ellipses.
                    Vc=orth([va vb]);
                    H=Vc'*A*Vc;
                    Y=nearritz(H,z,tol);
                    vf=Vc*Y;
                    rvf=vf'*A*vf;
                    break
                end
                
                %Set things up for the next iteration
                %Determine side of the triangle to which z is closest
                %Recall triangle points should be stored in
                %rva rvb and rvc
                %Sides of concern are rva rvc and rvb rvc
                
                cm=(rva+rvc+rvb)/3;
                %theangles=angle(([rva rvb rvc]-cm)*exp(-1i*angle(z-cm)));
                
                theangles=myunwrap(angle([rva rvc rvb]-cm));
                zangle= myunwrap([angle(rva-cm) angle(z-cm)]);
                zangle=zangle(2);
                theangles = theangles-zangle;
                
                an=sum(theangles<0);
                bn=an+1;
                if an==1
                    rvb=rvc;
                    vb=vc;
                else
                    rva=rvc;
                    va=vc;
                end
                jj=jj+1;
                if jj>iters
                    break
                end
                err0=err1;
                err1=abs(rvf-z);
                if abs(err0-err1)<tol
                    eeval=-eeval;
                    rvf=[];
                    vf=[];
                    break
                end
            end
        end
    else
        %    disp('Clearly not in the numerical range')
        rvf=[];
        vf=[];
        eeval=-4;
    end
end
warning on
return

    function [v]=eig_boundarypt
        [v s]=eig(HA);
        s=diag(s);
        [val I ]= max(s);
        v=v(:,I);
    end
    function [v]=eigs_boundarypt
        [v s]=eigs(HA,1,'LR',opts);
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ni fon]=insidequad(pc,p3,p4,tol)
%See Introduction to Algorithms Cormin Leiserson and Rivest
%Determines if the point p3 lies inside the polygon determined by pc.
%pc vertices of a  polygon in counterclockwise order.
%p3 point in question
%p4 point outside of polygon
%tol tolerance for how checking proximity

m=length(pc);
pc=[pc pc(1)];


x3=min(real(p3),real(p4));y3=min(imag(p3),imag(p4));
x4=max(real(p3),real(p4));y4=max(imag(p3),imag(p4));

ni=0;
fon=0;
for i=1:m
    p1=pc(1+(i-1));p2=pc(2+(i-1));
    x1=min(real(p1),real(p2));y1=min(imag(p1),imag(p2));
    x2=max(real(p1),real(p2));y2=max(imag(p1),imag(p2));
    if (x2-x3>=-tol) &&(x4-x1>=-tol)&&(y2-y3>=-tol)&&(y4-y1>=-tol)
        if ~((abs(p3-p1)<tol)||(abs(p3-p2)<tol))
            cp1=-imag((p3-p1)*conj(p2-p1))/abs(p3-p1)/abs(p2-p1);
            cp2=-imag((p4-p1)*conj(p2-p1))/abs(p4-p1)/abs(p2-p1);
            if (cp1*cp2<0)
                %The book failed to mention the possibility that
                %in this case p1p2 may not be long enough to intersect p3p4
                %Determine where they would intersect
                % and if this is on the segment p1p2
                
                angle1=abs(angle((p3-p1)*conj(p2-p1)));
                angle2=pi/2 ;
                
                x5=cos(angle1)*abs(p3-p1);
                z5=x5*exp(1i*angle(p2-p1))+p1;
                y5=sin(angle1)*abs(p3-p1);
                angle4=angle((p4-p3)*conj(z5-p3));
                z=(p4-p3)*sec(angle4)*y5/abs(p4-p3)+p3;
                
                
                if ~((abs(z-p1)<tol)||(abs(z-p2)<tol))%If p3p4
                    %passes
                    %through a corner
                    if (abs(angle((z-p1)*conj(p2-p1)))<pi/2)&&(abs(z-p1)<(1+1e-15)*(abs(p2-p1)))
                        
                        ni=ni+1;
                    end
                else
                    ni=ni+.5;
                    fon=fon+1;
                end
                
            else
                
                if ((abs(cp1)<1e-14 && cp2~=0) || (cp1~=0 && abs(cp2)<1e-14) ) || (((abs(cp1)<tol)&&(abs(cp2)<1e-14))&&(abs(p1-p2)>1e-15))
                    
                    %The p3 falls on the line p1p2 and our work is done
                    %fon=i;
                    ni=ni+1;
                end
                
            end
        else
            %p3 is a vertex of the polygon.
            ni=ni+.5;
            fon=fon+1;
        end
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%Need to find actual closest point
%Not closest on ray from center of ellipse.
%I have the semimajor and minor axis
%with everything oriented appropriatelly
%First we want to put the cz in the appropriate quadrant
cz2=(z-cH)*exp(-1i*rangle);

%%Need to adjust
%%Accounts for degenerate ellipse
if abs(b/a)<1e-14
    phi=acos(min( abs(cz2/d),1)*sign(real(cz2)/d));
else
    phi=atan2((imag(cz2)*d),(real(cz2)*b));
end
closetoz2 = d*cos(phi)+b*1i*sin(phi);

cz=closetoellipse(d,b,real(cz2),imag(cz2));
if abs(b/a)<1e-14
    phi=acos(min(abs(cz/d),1)*sign(real(cz)/d));
else
    phi=atan2((imag(cz)*d),(real(cz)*b));
end
closetoz = d*cos(phi)+b*1i*sin(phi);


if abs(cz2/closetoz2)<=1+tol
    % Point is in the ellipse
    if abs(b/a)<1e-14
        theta=pi/4;
    else
        phi=atan2((imag(cz2)*d),(real(cz2)*b));
        theta=  asin(cz2/closetoz2)/2;
    end
else
    theta=pi/4;
end
cz2=closetoz*exp(1i*rangle)+cH;
psi=1/2*atan(-a/b);
x=[sin(psi); cos(psi)];
x'*CH*x;
X=[sin(psi) cos(psi); cos(psi) -sin(psi)];
CHT = X'*CH*X;
if CHT(1,2)<0
    D=diag([1,-1]);
else
    D=eye(2);
end
CHTP=D'*CHT*D;
if CHTP(2,1)>CHTP(1,2)
    Dl=[0 1 ; 1 0];
else
    Dl=eye(2);
end
CHTP=Dl*CHTP*Dl;
CHTP=CHTP-diag(abs(diag(CHTP)));
y=[cos(theta); sin(theta)*exp(1i*phi)];
Y=U*Dm*X'*D'*Dl'*y;
pr=Y'*H*Y;
[ pr z];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cz=closetoellipse(a,b,u,v)
% Determines closest point on ellipse with semimajor and semiminor axis
%lengths of a and b centered at the origin to the point (u,v)
%See 'Distance from a Point to an Ellipse in 2D
%David Eberly http://www.geometrictools.com

%Modify (u,v), so it is in the first quadrant
quad=0;
if v<0
    quad=bitor(quad,2);
    v=-v;
end
if  u<0
    quad=bitor(quad,1);
    u=-u;
end

%Largest root of the following determines
%the closest point
tr= roots([(-1) (- 2*a^2 - 2*b^2) (a^2*u^2 - 4*a^2*b^2 - a^4 - b^4 + b^2*v^2) (2*a^2*b^2*u^2 - 2*a^2*b^4 - 2*a^4*b^2 + 2*a^2*b^2*v^2) (a^4*b^2*v^2 - a^4*b^4 + a^2*b^4*u^2)]);

t=tr(real(tr)>=-b^2*(1+1e-8));
[val ind]=sort(real(t),'descend');
%This wouldn't work in general
%however do to the nature of the inversefov algorithm
%the point will always have abs(u)<=a
if abs(v)<a*1e-14;
    t=0;
else
    t=val(1);
end
%f(a,b,u,v,t)
x=(a^2*u)/(t+a^2);
if b~=0
    y=(b^2*v)/(t+b^2);
else
    y=0;
end
if bitand(quad,2)
    y=-y;
end
if bitand(quad,1)
    x=-x;
end
cz=x+1i*y;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function x = myunwrap(y)
%This does what matlab's unwrap function should do
for i=1:(length(y)-1)
    while y(i+1)<y(i)
        y(i+1)=y(i+1)+2*pi;
    end
end
x=y;
end

