function b = izboljsava(A,Vbox,rvbox,tol)
%Ce je mozno ta funkcija najde izotropni vektor b, ki je resitev 0 = b'*A*b
%Vhodi: 	A 		matrika
%			Vbox	matrika najvecjih in najmansih lastnih vektorjev
%			rvbox	vektor robnih tock dobljenih kot x'*A*x, kjer so x vekotrji iz Vbox
%			tol		toleranca
%
%Izhod:		b		izotropni vektor
%
%Povzeto po R.Carden

q=real(rvbox(3))+ max(.1,.1*real(rvbox(3)-rvbox(1)))+1i*(imag(rvbox(4))+ max(.1,.1*imag((rvbox(4)-rvbox(2)))));

x1=min(real(rvbox));y1=min(imag(rvbox));
x2=max(real(rvbox));y2=max(imag(rvbox));

[ni fon] =insidequad([x1+1i*y1 x2+1i*y1 x2+1i*y2 x1+1i*y2],0,q,tol);

    if mod(ni,2)==1
        ni=0;
		
        [ni fon] =insidequad(rvbox,0,q,tol);
		
        if mod(ni,2)==1
		           angle1=angle((0-rvbox(1))*conj(rvbox(3)-rvbox(1)));
            
            if abs(angle1)<1e-14   %If 0 lies on rvbox(1)-rvbox(3)
                Vc=orth(Vbox(:,[1 3]));
                rvc=0;
            else
                if angle1>0 &&length(rvbox)~=3
                    rvbox=rvbox([1 3 4]);
                    Vbox=Vbox(:,[1 3 4]);
                end
                %Use law of sines to determine where ray from rvbox to 0 intersects.
                %ray from rvbox(3)to rvbox(2)
                angle1=angle((rvbox(3)-rvbox(1))*conj(0-rvbox(1)));
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
            
            %Determine a Ritz vector for 0
            if abs(rvc)<tol
                b=vc;
                rb=rvc;
            else
                
                Vc=orth([vc Vbox(:,1)]);
                H=Vc'*A*Vc;
                Y=nearritz(H,0,tol);
                b=Vc*Y;
                rb=b'*A*b;
            end
		end
	end
end
