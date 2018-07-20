function [D, dim, pars] = SNLExamples(problemname, noisetype, n, m, nf, range)
%
% problemname:  {'BLTWY06_Inner','Zh07','Tseng07','BNQ15','BLTWY06_Outer'} 
% collected from existing literature
%
% noisetype = 'additive', 'multiplicative', 'log-normal'
% n: number of anchors and sensors
% nf: 0.1
% range: used to construct the initial (perturbed) distances
%
% D = [anchor-anchor (squared) distance, anchor-sensor (squared) distance;
%      sensor-anchor (squared) distance, sensor-sensor (squared) distance]
%      distances are squared    
%      diag(D) = 0

rng('default');
randstate = rng('shuffle');


while 1
    
switch problemname        
    case 'BLTWY06_Inner'
        dim = 2;
        m   = 4; s = n - m;
        B   = 0.2*[1 1 -1 -1; 1 -1 1 -1]; 
        PA  = B;
        PS  = -0.5 + rand(2,s); 
        PP  = [PA, PS];
    case 'BLTWY06_Outer'
        dim = 2;
        m   = 4; s = n - m;
        B   = [0.45 0.45 -0.45 -0.45; 0.45 -0.45 0.45 -0.45];        
        PA  = B;
        PS  = -0.5 +rand(2,s);
        PP  = [PA, PS]; 
    case 'Zh07'
        dim = 2;
        m   = 5; s = n - m;
        B   = 0.4*[-1 1 1 -1 0; 1 1 -1 -1 0];
        PA  = B;
        PS  = -0.5 +rand(2,s); 
        PP  = [PA, PS];     
    case 'Tseng07'
        dim = 2;
        PP  = -0.5 + rand(2,n);
        PA  = PP(:, 1:m);
        PS  = PP(:, (m+1):end);                  
    case 'BNQ15'
        dim = 2;
        PP  = 100*rand(2,n);
        PA  = PP(:, 1:m);
        PS  = PP(:, (m+1):end);   
    case '3DSNL'
        dim = 3;
        PP  = rand(3,n);
        PA  = PP(:, 1:m);
        PS  = PP(:, (m+1):end);
    case '3DSNLsymetric'        
        dim = 3;
        m   = 8;
        PA  = [3 7 7 3 3 7 7 3;
               3 3 7 7 3 3 7 7;
               3 3 3 3 7 7 7 7]/10;
        PS  = rand(3,n-m);
        PP  = [PA PS];
    case 'Noanchors'
        dim = 2;
        m   = 0;
        PA  = [];
        PS  = 100*rand(2,n-m);
        PP  = [PA PS];     
    otherwise
        disp('input a problen name');
        
end


% %%%%% end of problemname

% generate the squared pre-distance matrix D
% use randidtance.m (Toh) to construct perturbed distances

D = randistance_qi(PA,PS,range,nf,noisetype,randstate);
% D = [DD, D0]
% DD -- sensor - sensor distances (pertrubed)
% D0 -- sensor - anchor distances (perturbed)

if m > 0
    D0 = squareform(pdist(PA'));
    D  = [D0 (D(:, (n-m+1):n))'; D(:, (n-m+1):n) D(1:(n-m), 1:(n-m))];
end

[~,bins] = graphconncomp(sparse(D));
if all(bins == 1);   break;
else  error('Graph is not connected, try to increase range...')
end

end

pars.m       = m;
pars.range   = range;
pars.PP      = PP;

end



% This is the code from Toh's SNLSDP package
% The only difference is that the multiplicative perturbation 
% used abs 9used in Tseng's paper (Tseng SIOPM 2007).
%
%%*************************************************************************
%% Compute pair-wise distances within R limit between known-unknown and
%% unknown-unknown pairs.
%%
%% Dall = randistance(P0,PP,Radius,nf,noisetype);
%%
%% P0 : anchor positions
%%      nfix being the number of anchors
%% PP : sensor positions 
%%      npts being the number of unknown sensors
%% nf : noise factor 
%% noisetype : 1 - Normal :  dnoisy = dactual + N(0,1)*nf (Default)
%%             2 - Multiplicative Normal : dnoisy= dactual*(1 + N(0,1)*nf)
%%             3 - Log Normal : dnoisy= dactual * 10^(N(0,1)*nf)
%% Dall = [DD, D0]
%% DD = [sensor-sensor distance];
%% D0 = [anchor-sensor distance];
%%*************************************************************************

  function [Dall] = ...
            randistance_qi(P0,PP,Radius,nf,noisetype,randstate)

  if ~exist('noisetype'); noisetype = 'additive'; end
  if ~exist('randstate'); randstate = 0;
      rng('shuffle');
  end

%  randn('state',randstate);
%  rand('state',randstate);
%  rng(randstate);
  
  [dim,npts] = size(PP);
  nfix = size(P0,2);
  D0 = sparse(npts,nfix);
  DD = sparse(npts,npts); 
%%
  if strcmp(noisetype,'additive'); noisetype = 1; end
  if strcmp(noisetype,'multiplicative'); noisetype = 2; end
  if strcmp(noisetype,'log-normal'); noisetype = 3; end
%%
  for j = 1:npts    
     if (nfix > 0)
        tmp = PP(:,j)*ones(1,nfix)-P0;
        rr = sqrt(sum(tmp.*tmp));      
        idx = find(rr < Radius);
        rr = rr(idx); 
        if (~isempty(idx))
            if (noisetype == 1)
                rr = rr + (randn(1,length(idx))*nf);
            elseif (noisetype == 2)
                rr = rr .*((1+(randn(1,length(idx))*nf))); % abs used
            elseif (noisetype == 3)
                rr = rr .* 10.^(1+(randn(1,length(idx))*nf));
            end
            D0(j,idx) = rr';
        end
     end        
     if (j > 1)
        tmp = PP(:,j)*ones(1,j-1) - PP(:,1:j-1);
        rr = sqrt(sum(tmp.*tmp));
        idx = find(rr < Radius);
        rr = rr(idx); 
        if (~isempty(idx))
            if (noisetype == 1)
                rr = rr + (randn(1,length(idx))*nf);
            elseif (noisetype == 2)
                rr = rr .*((1+(randn(1,length(idx))*nf))); % abs used
            elseif (noisetype == 3)
                rr = rr .* 10.^(1+(randn(1,length(idx))*nf));
            end
            DD(idx,j) = rr';
        end
     end
  end
  DD = triu(DD,1) + triu(DD,1)';
  Dall = [DD, D0];
  end
  
 