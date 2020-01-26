function Out = PREEEDM(D,dim,pars)
%
% This code aims to solve the model
%
%   min_Z  \| H.*(sqrt(Z)-D) \|_1 + (rho/2) \|Z+P_Kr(-Z)\|^2
%    s.t.    L<=Z<=U
%
%
% INPUTS:
%
%   D   : the (n x n) dissimilarities matrix                     [required]
%         diag(D) = 0
%         dissimilarities are UNSQUARED, i.e.,
%                          D_ij=||point_i-point_j||_2+noise 
%         
%   dim : the embedding dimension  (e.g. dim = 2 or 3)           [required]
%
%   pars: parameters and other information                       [optional]
%         pars.m   : m  -- the number of given points, m>=0  
%         pars.PP  : PP -- (dim x n) matrix of coordinates of n points with
%                           first m(>=0) columns being given
%                    For sensor network localization (SNL)
%                    PP = [PA, PS]
%                    PA -- (dim x m) matrix of coordinates of m anchors
%                    PS -- (dim x (n-m)) matrix of coordinates of (n-m) sensors
%        pars.rho  : penalty parameter, default one rho = sqrt(n)*nnz(D)*max(D(:))/n^2
%        pars.LOWBD: lower bound i.e., L=pars.LOWBD.^2, Z_{ij}>=L_{ij}
%        pars.UPPBD: upper bound i.e., U=pars.LOWBD.^2, Z_{ij}<=U_{ij}
%                    Note: elements of pars.LOWBD and pars.UPPBD are UNSQUARED distances                          
%        pars.range: the communication range of two points, which means
%                    upper bound for Z_{ij}<=pars.range^2 if (D_{ij}>0  & i~=j)
%                    lower bound for Z_{ij}>=pars.range^2 if (D_{ij}==0 & i~=j)
%                    Note: pars.range is particular for SNL problem. If pars.range
%                    exists, no need pars.LOWBD and pars.UPPBD
%        pars.Otol : tolerance for objective,  default Otol=nnz(D)*1e-4 
%        pars.Etol : tolerance for eigenvalue, default Etol=1e-2 
%        pars.draw : 1--plot localizations in Re^dim  
%                    0--no plot (default) 
%
% OUTPUTS:
%
% If NO pars.PP exists 
%       Out.X:        dim-by-n matrix,  final coordinates 
%       Out.Time:     total time
%       Out.f:        relative f/f(0)
%
% If pars.PP exists 
%       Out.X:        (dim x (n-m)) matrix, coordinates before refinement 
%       Out.rX:       (dim x (n-m)) matrix, coordinates after refinement 
%       Out.Time:     total time including time for refinement
%       Out.rTime:    time for refinement
%       Out.RMSD:     Root Mean Square Distance (RMSD) before refinement 
%       Out.rRMSD:    Root Mean Square Distance (RMSD) after refinement 
%
% Refinement step is taken from Kim-Chaun Toh SNLSDP solver
%
% Send your comments and suggestions to   [ sz3g14@soton.ac.uk ]                       
%
% Warning: Accuracy may not be guaranteed!!!!!   
%
% This version: March 1st, 2018,   written by Shenglong Zhou    

t0=tic;
% parameters design -------------------------------------------------------
n     = size(D,1);
nD    = nnz(D);
rate  = nD/n/n;
if nargin==2; pars=[]; end
[m,itmax,Eigtol,Objtol] = getparameters(nD,pars);
if m>0  
    fprintf('\nNumber of given points : %3d\n',m);
    fprintf('Number of unknown points: %3d\n',n-m);
    fprintf('Procrustes analysis and refinements step will be done!\n');
else
    fprintf('\nNo points are given!\n'); 
    fprintf('Number of unknown points: %3d\n',n);
end



% shortest path to complete the missing elements of D----------------------
    fprintf('Available dissimilarities rate: %1.2f \n',rate);
if  rate<0.05
    fprintf('Suggest providing more available dissimilarities!\n');
end
if rate<0.9
    DD = D.*D;        
    fprintf('\nContruct the shortest paths...\n');
    ts=tic;
    SD = graphallshortestpaths(sparse(sqrt(DD)));
    if any(SD == Inf)
    error('The neighborhood graph is not connected, increase range!');
    end
    SD = max(SD, D);  % to only replace those missing distances
    fprintf('Done the shortest paths by using %1.2f seconds\n',toc(ts));
else
    SD = D;
   fprintf('No shortest paths calculated!\n');
end

% scale the elements of D ---------------------------------------------
fSD   = full(SD);
scale = max(fSD(:)); 
Do    = D;
if    scale<=10; scale=1; 
else; D=D./scale; fSD=fSD./scale; 
end

D  = full(D);
DD = D.*D; 
H  = full(spones(D));
r  = dim;
T  = [];
if m>0; T = 1:m; H(T,T)=0;  end

% ilitialize ( L  U  Z c ) ------------------------------------------------
Z  = fSD.*fSD; 
UB = n*max(Z(:));
L  = zeros(n); 
U  = UB*ones(n);

if isfield(pars,'LOWBD')   
   L  = (pars.LOWBD/scale).^2;  
end

if isfield(pars,'UPPBD')   
   U  = (pars.UPPBD/scale).^2;
   HU = spones(U);
   if nnz(HU)<(n^2-n); U=U+(1-HU).*UB; end
   if max(U(:))==inf;  U(U==inf)=UB;   end 
end 

if isfield(pars,'range')   
   H1 = 1-H; 
   rs = (pars.range/scale)^2;
   L  = L.*H  + H1*rs;                           
   U  = U.*H1 + H*rs;
end

L(T,T)       = Z(T,T);      
U(T,T)       = Z(T,T); 
L(1:n+1:end) = 0;               
U(1:n+1:end) = 0;  
 
Z      = min( U,  max(L, Z ));  
rho    = sqrt(n)*rate*max(D(:));%sqrt(nD)/log(n)/max(D(:))
if isfield(pars,'rho')  
rho    = pars.rho; 
end
 
Hr     = H/rho;
TH     = find(Hr>0);  
PZ     = ProjKr(-Z,r);
frZ    = sum(sum(abs(sqrt(Z(TH))-D(TH)))) +(rho/2)*FNorm(Z+PZ);
OErr   = zeros(itmax,1);
EErr   = zeros(itmax,1);
fprintf('Start to run ... \n');
fprintf('\n------------------------------------------------\n');

% main loops --------------------------------------------------------------
for iter=1:itmax
  
  
    Z   = dcroot( -PZ, Hr, TH, DD,L, U );       
    PZ  = ProjKr(-Z,r);
    
    % stop criteria
     gZ      = FNorm(Z+PZ);
     ErrEig  = gZ/FNorm(JXJ(-Z));
     frZo    = frZ;      
     frZ     = sum(sum(abs(sqrt(Z(TH))-D(TH)))) +(rho/2)*gZ;
     ErrObj  = abs(frZo-frZ)/(1+rho+frZo);
     fprintf('Iter: %3d  ErrEig: %.3e  ErrObj: %.3e\n',iter, ErrEig, ErrObj);
 
     if  iter>=5 && ErrEig<Eigtol && ErrObj<Objtol; break; end
     
     % update rho
     OErr(iter,1)=ErrObj; EErr(iter,1)=ErrEig;
     if  ErrEig>Eigtol && ErrObj < Objtol/5&&...
         iter>=10 && var(EErr(iter-9:iter,1))/ErrEig<1e-4
         rho = 1.25*rho; Hr = H/rho; 
         Objtol = min(Objtol, max(Objtol/1.1, ErrObj));
     end
     if  ErrObj>Objtol && ErrEig<Eigtol/5  
         if iter<=5
         rho = 0.5*rho; Hr = H/rho;
         elseif iter>5 && var(OErr(iter-4:iter,1))/ErrObj<1e-4
         rho = 0.75*rho; Hr = H/rho;
         end
     end
 
      
end

% Results and graph output -----------------------------------------------
if isfield(pars, 'PP')  
    T1            = 1+m:n;
    PPs           = pars.PP/scale;
    Xs            = procrustes_zhou(m, PPs, Z);  
    Out.Time      = toc(t0); ts=tic;   
    rX            = refinepositions(Xs,PPs(:,T),[D(T1,T1) D(T1,T)]);%refinement step
    Out.rTime     = toc(ts);  
    Out.Z         = Z*scale^2;
    Out.X         = Xs*scale;
    Out.rX        = rX*scale;
    Out.Time      = Out.Time+Out.rTime;
    Out.RMSD      = sqrt(FNorm(pars.PP(:,T1)-Out.X)/(n-m));
    Out.rRMSD     = sqrt(FNorm(pars.PP(:,T1)-Out.rX)/(n-m));   
    if isfield(pars, 'draw') && pars.draw
        figure
        subplot(211)
        plot_SNL(Out.X,Out.RMSD,0,pars);  
        subplot(212)
        plot_SNL(Out.rX,Out.rRMSD,1,pars); 
    end
    fprintf('------------------------------------------------\n');
    fprintf('rTime:     %1.3fsec\n',  Out.rTime)
    fprintf('Time:      %1.3fsec\n',  Out.Time)
    fprintf('RMSD:      %1.2e \n',    Out.RMSD )
    fprintf('rRMSD:     %1.2e \n',    Out.rRMSD)
else
    Out.Time      = toc(t0);
    [U,E]         = eig(JXJ(-Z)/2);
    Eig           = sort(real(diag(E)),'descend')*(scale^2);
    Er            = real((Eig(1:r)).^0.5); 
    Out.Eigs      = Eig;
    Out.X         = sparse(diag(Er))*real(U(:,1:r))';
    Out.f         = sqrt(FNorm(sqrt(Z(H~=0))-Do(H~=0))/FNorm(Do(H~=0)));
    fprintf('------------------------------------------------\n');
    fprintf('Time:      %1.3fsec\n',     Out.Time)
    fprintf('Stress:    %1.2e \n',    Out.f)    
end

end

% ------------------------------------------------------------------------
function  [m,itmax,Etol,Otol] = getparameters(nD,pars)
itmax  = 2000; 
m      = 0;
Otol   = log(nD)*1e-4;
Etol   = 1e-2;
if isfield(pars, 'Otol');   Otol=pars.Otol;   end
if isfield(pars, 'Etol');   Etol=pars.Etol;   end
if isfield(pars, 'm');      m   = pars.m;     end

end
 

% ------------------------------------------------------------------------
function  JXJ = JXJ( X )
% Compute J*X*J where J = I-ee'/n;
    nX   = size(X,1);
    Xe   = sum(X, 2);  
    eXe  = sum(Xe);    
    JXJ  = repmat(Xe,1,nX);
    JXJt = repmat(Xe',nX,1);
    JXJ  = -(JXJ + JXJt)/nX;
    JXJ  = JXJ + X + eXe/nX^2;
end


% -------------------------------------------------------------------------
function fn= FNorm(A)
% Compute the Frobenius norm of A, i.e., ||A||_F^2
    fn=sum(sum(A.*A));
end

% -------------------------------------------------------------------------
function Z0= ProjKr(A,r)
% The projection of A on cone K_+^n(r)  
        JAJ    = JXJ(A);
        JAJ    = (JAJ+JAJ')/2;
        [V0,P0]= eigs(JAJ,r,'LA');
        Z0     = real(V0*max(0,P0)*V0'+A-JAJ); 
end

% -------------------------------------------------------------------------
function Xs= procrustes_zhou(m0, PP, D0)
% Input:
% m0 --the number of given anchors in Re^r0
% PP --r0-by-n0 matrix whose first m0 columns are  given anchors 
%      and rest (n0-m0) columns are given sensors
% D0 --n0-by-n0 distance matrix 
%     [anchors-anchors (squred) distance,anchors-sensors  (squred) distance
%      sensors-anchors (squred) distance,sensors-sensors  (squred) distance]
%      containing m0 computed anchors and (n0-m0) sensors in Re^r0
% Output:
% Xs -- r0-by-(n0-m0) matrix containing (n0-m0) sensors in Re^r0
    r0     = size(PP,1);
    n0     = size(D0,2);
    JDJ    = JXJ(D0);
    JDJ    = -(JDJ+JDJ')/4;
    [U1,D1]= eigs(JDJ,r0);   
    X0     = (D1.^(1/2))*(U1');   
    if m0>0
    A      = PP(:,1:m0);
   [Q, ~, a0, p0] = procrustes_qi(A,X0(:,1:m0));	
    Z0     = Q'*(X0-p0(:, ones(n0, 1))) + a0(:, ones(n0, 1)); 
    Xs     = Z0(:,m0+1:n0);
    Xa     = Z0(:, 1:m0);
    %Xs    = Xs*mean(sum(Xa.*A)./sum(Xa.^2));
    Xs     = Xs*max(1,sum(sum(Xa.*A))/sum(sum(Xa.^2)));
    else        
    [~,~,Map]= procrustes(PP',X0'); 
    Xs       = (Map.b*Map.T')*X0 + diag(Map.c(1,:))*ones(size(X0));
    end
end

% -------------------------------------------------------------------------
function [Q, P, a0, p0] = procrustes_qi(A, P)
% Procrustes analysis for rigid localization of anchors in A
% Input:
% A -- r-by-m0 matrix containing m0 anchors in Re^r
% P -- r-by-m0 matrix containing m0 computed locations of m0 anchors in Re^r
% 
% Output:
% Q  -- the orthogonal transformation
% P  -- the resulting coordinators of P after transformation
% a0 -- the tranlation vector that simply centrizes A
% p0 -- the tranlation vector that simply centrizes P
%
%
 m0 = size(A,2);
% centerizing A and P
 a0 = sum(A, 2)/m0;
 p0 = sum(P,2)/m0;
 A  = A - a0(:, ones(m0,1));
 P1 = P - p0(:, ones(m0,1));
 P  = P1*A';
 [U0, ~, V] = svd(P);
 Q  = U0*V';
 P  = Q'*P1 + a0(:, ones(m0,1));
end
