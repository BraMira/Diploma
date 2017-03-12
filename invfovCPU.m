function [x err step] = invfovCPU(A,u,gr,pr) % copyright Frank Uhlig, Feb. 2010, revised March 2011
% Input :  A is an n by n matrix;
%          u is a given point in the complex plane for which we seek a unit vector x with x'*A*x = u 
%          if gr = 1 we plot the approximate boundary of the FOV of A, as well as computed points
%          if pr = 1 we give intermediate output on-screen
% Output : x is a generating complex unit vector for u in the field of values of A, i.e., x'*A*x = u 
%          err = norm(x'*A*x - u), almost zero
%          step counts the eigen evaluations that were needed to find x (if possible)
%          If u lies outside FOV(A) then  x = zero vector and  err  = inf
% To accompany the paper:
%    Ch. Chorianopoulos, P. Psarrakos, and F. Uhlig, A method for the inverse numerical range problem,
%    Electronic Journal of Linear Algebra (ELA), (20) 2010, p. 198 - 206.
% Revised in March 2011 thanks to comments by Gerard Meurant (index mess-up on line 98) and 
%    by  incorporating the stable version of the 'quadratic formula' in all
%    instances (lines 110, 122)

if nargin == 1, u = 0; gr = 0; pr = 0; end           % preliminaries:
if nargin == 2, gr = 0; pr =  0; end
if nargin == 3, pr = 0; end
[n,m] = size(A); if n ~= m, frintf(' Matrix A is not square: we quit!'); return, end
A = A-u*eye(n);                                      % shift problem to zero in F(A-uI)
x = zeros(n,1); z = x; err = inf; maxit = 40; d = 0; % start-up values
opts.disp = 0; opts.maxit = 1000; opts.tol = 10^-4; warning('OFF'); kryln = 200; rpd = []; rmd = [];
H = (A+A')/2; K = (A-A')/(2i);                       % Hermiten/skewhermitean parts H, K of A = H + iK
if gr == 1, wber3(A,60,opts,gr,kryln); hold on, end  % plot FOV boundary if desired

%    Find four boundary FOV points rm, rM on left/right and im, iM on bottom, top and their generators
[rM,rm,x1p,x1n,step,x,err,d] = starteig(A,H,n,opts,kryln,u,gr,x,err); % two FOV bdry pts rm, rM
if d == 1, if gr == 1, hold off, end, return, end                     % done if zero is encircled
[iM,im,y1p,y1n,s,x,err,d] = starteig(A,K,n,opts,kryln,u,gr,x,err); step = step + s; % FOV bdry pts im,iM
if d == 1, if gr == 1, hold off, end, return, end                     % done
%    Find real axis points of great circle images through any 2 points among im, iM, rm, rM
[rpd,rmd,x,err,d] = realposnegdat(A,iM,y1p,im,y1n,u,gr,rpd,rmd);      % iM-im segment real ellipse pts
if d == 1, if gr == 1, hold off, end, return, end                     % done
if rm(2) < 0 
  [rpd,rmd,x,err,d] = realposnegdat(A,iM,y1p,rm,x1n,u,gr,rpd,rmd);             % ditto for iM-rm segment
  if d == 1, if gr == 1, hold off, end, return, end                            % done
else 
  [rpd,rmd,x,err,d] = realposnegdat(A,im,y1n,rm,x1n,u,gr,rpd,rmd);             % or for im-rm segment
  if d == 1, if gr == 1, hold off, end, return, end, end                       % done
if rM(2) < 0, [rpd,rmd,x,err,d] = realposnegdat(A,iM,y1p,rM,x1p,u,gr,rpd,rmd); % ditto for iM-rM segment
  if d == 1, if gr == 1, hold off, end, return, end                            % done
else [rpd,rmd,x,err,d] = realposnegdat(A,im,y1n,rM,x1p,u,gr,rpd,rmd);          % or for im-rM segment
  if d == 1, if gr == 1, hold off, end, return, end, end                       % done
if rm(2)*rM(2) < 0, [rpd,rmd,x,err,d] = realposnegdat(A,rm,x1n,rM,x1p,u,gr,rpd,rmd); % last: rm-rM 
  if d == 1, if gr == 1, hold off, end, return, end, end                             % done
if isempty(rpd) == 0 && isempty(rmd) == 0   % if FOV bdry points lie on both sides of zero on real axis:
  [~, k] = max(rmd(1,:)); [~, j] = min(rpd(1,:)); k = k(1); j = j(1);  
  [x err] = realaxisintpol(A,rmd(1,k),rpd(1,j),rmd(2:end,k),rpd(2:end,j),u);   % done
  if gr == 1, hold off, end, return, end

%    If not done, find closest to zero FOV point on one side of zero and set up for iteration
if isempty(rpd)     % set up closest real axis point for iteration and their sustaining angles
  [~, k] = max(rmd(1,:)); k = k(1); R = rmd(1,k); vR = rmd(2:end,k); end % only neg. real axis pts
if isempty(rmd)     % [~,k] = min(..) works only from MATLAB 2009b on; earlier versions use [nn,k] = ..
  [~, k] = min(rpd(1,:)); k = k(1); R = rpd(1,k); vR = rpd(2:end,k); end % only pos. real axis pts
if R(1) < 0, P = rM; vP = x1p; h1 = [1,0]; else P = rm; vP = x1n; h1 = [-1,0]; end % sustaining 
if P(2) < 0, Q = iM; vQ = y1p; k1 = [0,1]; else Q = im; vQ = y1n; k1 = [0,-1]; end %    angles
ST = step; st = 0;                         % step counts eigen evaluations
%    Iterate by bisection of the sustaining angle until a decision can be made
while st < maxit    
  ra = (h1 + k1)/2; ra = ra/norm(ra); At = ra(1)*H + ra(2)*K;   % bisect angle
  vN = eigsw(At,n,opts,kryln,u); st = st + 1;                   % compute new FOV bdry point N
  if pr == 1, fprintf('iterating: at step  st = %g  \n',st), end
  if norm(vN) == 0, x = z; err = inf; step = ST + st; if gr == 1, hold off, end, return, end % u not in FOV
  N = vN'*A*vN; N = [real(N);imag(N)];     % new FOV point N with generating vector vN
  if gr == 1, plot(N(1),N(2),'dr'),  fprintf('Iter \n'), pause, end
  if N(2)*P(2) < 0                         % if N lies opposite P find real ellipse points through N and P
    X = Cvectinterpol(A,P,vP,N,vN,u); xAx = real(diag(X'*A*X)); % find real axis points
    if gr == 1, plot(real(xAx),zeros(length(xAx),1),'+k'), end
    j = find(R(1)*real(xAx) < 0); 
    if isempty(j), k1 = ra; Q = N; vQ = vN;                     % update boundary points           
    else                                                        % find generator x for u = x'*A*x
      j = j(1); [x err] = realaxisintpol(A,R,xAx(j),vR,X(:,j),u); step = ST + st; 
      if gr == 1, hold off, end, return, end
  else P = N; vP = vN; h1 = ra; end, end   % or just update sustaining FOV bdry pts

 %   Seven subroutines used in  invfovCPU.m :   
 function [Ma,Mi,xMa,xMi,step,x,err,d] = starteig(A,H,n,opts,kryln,u,gr,x,err) % two FOV bdry points
  [xMa xMi s d] = eigsw(H,n,opts,kryln,u); step = s; d = 0;     % extreme eigenvectors of H
  if norm(xMa) == 0,                                            % u not in FOV(A): stop
    if gr ==1, hold off, end, err = inf; x = xMa; d = 1; Ma = inf; Mi = Ma; return, end 
  Ma = xMa'*A*xMa; Ma = [real(Ma);imag(Ma)];                    % evaluate xMa'*A*xMa 
  if gr == 1, plot(Ma(1),Ma(2),'sr'), end        
  if norm(Ma) < 10^-13, err = norm(Ma); x = xMa;                % approximate generator for u found:
    if gr == 1, hold off, end, d = 1; return, end               % Done !
  Mi = xMi'*A*xMi; Mi = [real(Mi);imag(Mi)];                    % evaluate x'*A*x 
  if gr == 1, plot(Mi(1),Mi(2),'sb'), end         
  if norm(Mi) < 10^-13, err = norm(Mi); x = xMi; d = 1;         % approximate generator for u found:
    if gr == 1, hold off, end, return, end                      % Done !

 function [rp,rm,x,err,d] = realposnegdat(A,P,vP,Q,vQ,u,gr,rp,rm) % real axis pts on gr. circle ellipse
   X = Cvectinterpol(A,P,vP,Q,vQ,u); xAx = real(diag(X'*A*X));    % call great circle interpolator   
   if gr == 1, plot(real(xAx),zeros(length(xAx),1),'+k'), end
   k = find(xAx > 0); l = find(xAx < 0);                        % separate pts into right/left side sets
   if isempty(k) == 0 && isempty(l) == 0, d = 1;                % if points on both sides of zero
     r1 = xAx(l(1)); r2 = xAx(k(1));   % indices for X.. below corrected by Gerard Meurant in March 2011
     [x err] = realaxisintpol(A,r1,r2,X(:,l(1)),X(:,k(1)),u);   % Done!  
     if gr == 1, hold off, end, return,   
   else d = 0; x = zeros(length(A),1); err = inf;               % else collect in right/left of zero sets
     rp = [rp, [xAx(k)';X(:,k)]]; rm = [rm, [xAx(l)';X(:,l)]]; end %   and add to data sets rm and rp
    
function [x err] = realaxisintpol(A,aa,cc,xaa,xcc,u)            % find a generator x of 0 in F(A)
  phi = angle(xcc'*A*xaa - xaa.'*conj(A)*conj(xcc));            %    ala Horn + Johnson, Topics, p.25
  alphminusphi = real(exp(1i*phi)*xaa'*A*xcc + exp(-1i*phi)*xcc'*A*xaa);
  discr = alphminusphi^2 - 4*aa*cc;
  if discr < 0,
    fprintf('[1]  u = %.18g%+.18gi  lies outside the FOV(A) (real interp *)\n',real(u),imag(u)), 
    x = []; err = inf; return, end
  t1 = real((-alphminusphi - mysign(alphminusphi)*discr^.5)/(2*cc)); t2 = aa/(cc*t1);
  x = exp(-1i*phi)*xaa + t1*xcc; x = x/norm(x); err = norm(x'*A*x); % two solutions:
  x2 = exp(-1i*phi)*xaa + t2*xcc; x2 = x2/norm(x2); err2 = norm(x2'*A*x2);
  if err2 < err, x = x2; err = err2; end                        % take smallest error x

function X = Cvectinterpol(A,Q,vQ,P,vP,u)               % find generators on great circle through vQ, vP
  R = vQ'*A*vP + vP'*A*vQ; R = [real(R),imag(R)];       %        that create FOV points on the real axis
  q = Q(2); r = R(2); p = P(2); f = p+q-r; pof = p/f;   % ala Ch. Davis, 1971 Can Math Bull
  go2 = r/(2*f) - pof;
  if go2^2 - pof < 0,                                   % first via tx + (1-t)y from x to y :
    fprintf('[2]  u = %.18g%+.18gi  lies outside the FOV(A) (complex interp *)\n',real(u),imag(u)) 
    X = []; return, end
  t1 = -go2 - mysign(go2)*(go2^2 - pof)^.5; t2 = pof/t1;% stable use of quadratic formula
  x1 = t1*vQ + (1-t1)*vP; x1 = x1/norm(x1); 
  x2 = t2*vQ + (1-t2)*vP; x2 = x2/norm(x2); X = [x1,x2];
   
function [vR vS step d] = eigsw(Atest,n,opts,kryln,u)   % extreme eigenvetors via QR or Krylov 
  if nargout == 1, step = 1; d = 0;                     % switch between start-up and iterations 
    if n >= kryln, [vR,ev] = eigs(Atest,1,'lr',opts);   % Use Krylov method for n >= kryln
      if ev < 0                                         % done
        fprintf('[3] One of H or K is definite and  u = %.18g%+.18gi  lies outside FOV(H + iK)\n',...
                 real(u),imag(u)), vR = zeros(n,1); d = 1; return, end       
    else [V,ev] = eig(Atest);                           % otherwise (n < kryln) use QR algorithm :
      if ev(1)*ev(end) > 0                              % done
        fprintf('[3] One of H or K is definite and  u = %.18g%+.18gi  lies outside FOV(H + iK)\n',...
                 real(u),imag(u)), vR = zeros(n,1); d = 1; return
      else vR = V(:,n); end, end 
  else d = 0;                                           % used at start-up
    if n >= kryln, step = 2;                            % Use Krylov method for n >= kryln
      [vR,ev1] = eigs(Atest,1,'lr',opts); [vS,ev2] = eigs(Atest,1,'sr',opts); % use eigs  
      if ev1*ev2 > 0                                    % done
        fprintf('[4] One of H or K is definite and  u = %.18g%+.18gi  lies outside FOV(H + iK)\n',...
                 real(u),imag(u)), vR = zeros(n,1); vS = vR; d = 1; return, end
    else step = 1;                                      % otherwise use QR algorithm :
      [V,ev] = eig(Atest); 
      if ev(1)*ev(end) > 0                              % done
        fprintf('[4] One of H or K is definite and  u = %.18g%+.18gi  lies outside FOV(H + iK)\n',...
                 real(u),imag(u)), vR = zeros(n,1); vS = vR; d = 1; return
      else vR = V(:,n); vS = V(:,1); end, end, end 
  
function wber3(A,L,opts,gr,kryln)                       % plot field of values boundary :
  n = length(A); PP = NaN * ones(2*(L+1),1);            % definitions 
  H = (A+A')/2; K = (A-A')/(2i);                        % hermitean and skewherm parts
  th = 0:pi/L:pi; ko = cos(th); si = sin(th);           % setting up angular values
  for k = 1:L+1,                                        % rotated evalue computations
    Ath = ko(k) * H + si(k) * K; 
    if n < kryln                                        % for small n < kryln use QR algorithm
      [V D] = eig(Ath);    
      Vt = V(:,n); PP(k) = Vt'*A*Vt;                    % desired evalues/evectors
      Vt = V(:,1); PP(L+1+k) = Vt'*A*Vt;                % (use max and min evalue/vector)
    else                                                % for large n use Krylov methods
      [V D] = eigs(Ath,1,'Lr',opts); PP(k) = V'*A*V;    % use eigs 
      [V D] = eigs(Ath,1,'Sr',opts); PP(L+1+k) = V'*A*V; end, end
  if gr == 1, plot(real(PP),imag(PP),'k'); hold on, grid on,    % set 1:1 aspect ratio for plot 
    plot(0,0,'or'), set(gca,'PlotBoxAspectRatio',[1,1,1]);  set(gca,'DataAspectRatio',[1,1,1]);  
    title(['Numerical range of a ',num2str(n),' by ',num2str(n),' matrix A'],'FontSize',14), 
    xlabel('real  axis', 'FontSize',12), ylabel('imaginary  axis', 'FontSize',12), end, 

function s = mysign(x)                        % sign function for real x that does not return zero
  s = 1; if x < 0, s = -1; end