function fovals(A,k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Author: Wilmer Henao    wi-henao@uniandes.edu.co
%   Department of Mathematics
%   Universidad de los Andes
%   Colombia
%
%   The Field of values of a matrix is a convex set in the complex plane
%   that contains all eigenvalues of the given matrix, this m-file plots 
%   boundary points of this set and the eigenvalues of the given matrix.
%
%   fovals(A,k)
%   A = The square matrix
%   k = The number of steps (500 by default) you can call FoV, without
%   this argument
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

EPS = 0.0000000000000000000001*i;
[m,n]=size(A);
if m ~= n
   disp('The matrix must be square')
   return
end

[Vo,Do] = eig(A);

tr = trace(A)/m;
A = A - tr.*eye(m);

if nargin == 1
    k = 500;
end

for ind = 1:k
    theta = 2*pi*ind/k;
    [vect,D] = eig(0.5.*((exp(i*theta).*A)+(exp(i*theta).*A)'));
    [a,b] = max(D*ones(m,1));
    V = vect(:,b)./norm(vect(:,b));
    pto(ind) = V'*A*V;
end

plot(pto + (tr*ones(1,k)) + (EPS*ones(1,k)),'k');
hold on
plot(Do*ones(m,1) + (EPS*ones(m,1)),'x','MarkerSize',6, 'MarkerEdgeColor', 'k', 'MarkerFaceColor','g', 'LineWidth',1)
% hold on
% plot(1,3,'o')
% hold on
% plot([-0.5910,-0.5910,2.9460,2.9460],[-3.0886,3.0886,- 3.0886,3.0886],'sk')
% hold on
% plot([-0.5910;-0.5910;2.9460],[-3.0886;3.0886;3.0886],'-k')
% hold on
% plot([2.9460;2.9460;-0.5910],[3.0886;-3.0886;-3.0886],'-k')
% hold on
% plot([2.9460, -0.5910,0.5415,0.5415],[0,0,3.0886,-3.0886],'dm','MarkerSize',10)
% hold on
% plot([2.9460; 0.5415;-0.5910;0.5415;2.9460],[0;3.0886;0;-3.0886;0],'--m')
hold off




% B=gallery('grcar',32);
% eeval=0; %Number of eigenvalue evaluations performed
% Vbox=zeros(32,4);
% rvbox=zeros(1,4);
% 
% %Determine four points on the boundary of numerical range
% for i=1:2
%     th=((i-1)*pi/2)-pi; %kako tece theta %od -pi=0? do -pi/2?
%     AA=exp(-1i*th)*(B); %e na i*theta *A %
%     HA=(AA+AA')/2; %H_theta
%     [va sa]=eig(HA);eeval=eeval+1; %va=l.vekt, sa=l.vrednosti
%     sa=diag(sa); %zloži l.vred v en stolpec
%     [val I ]= max(sa); %val najveèja l. vred, I na katerem mestu je 
%     Vbox(:,1+(i-1))=va(:,I); %v stolpec da na I to mesto, I to l. vred
%     rvbox(1+(i-1))=va(:,I)'*B*va(:,I); %dela vektor l. vred pomn z A
%     [val I ]= min(sa);
%     Vbox(:,3+(i-1))=va(:,I);
%     rvbox(3+(i-1))=va(:,I)'*B*va(:,I);
% end
% 
% 
% %pick a point outside the numerical range for determining if $z$ is in
% %polygonal approximation
% q=real(rvbox(3))+ max(.1,.1*real(rvbox(3)-rvbox(1)))+1i*(imag(rvbox(4))+ max(.1,.1*imag((rvbox(4)-rvbox(2)))));
% 
% %Form outer approximation to numerical range
% x1=min(real(rvbox));y1=min(imag(rvbox));
% x2=max(real(rvbox));y2=max(imag(rvbox));
% obox=[x1+1i*y1 x2+1i*y1 x2+1i*y2 x1+1i*y2]


% [x err step] = invfovCPU(B,0,1,1) z odkomentiranjem dveh ukazov