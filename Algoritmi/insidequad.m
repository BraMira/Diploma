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
