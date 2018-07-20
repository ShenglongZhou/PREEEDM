clc; clear; close all;  
 
%Initialization% ----------------------------------------------------------
test    = 1;
proname = {'1LFB.mat','1RGS.mat'};
load(proname{1})  

Inputs.noiseType     = 'normal';%uniform 
Inputs.sparseLevel   = 0.5;     % Note: use 0.3-1.0 to make problem easier
Inputs.noiseLevel    = 0.1; 
Inputs.radius        = 6;
Inputs.nDimensions   = 3;
Inputs.minLowerBound = 1.0;     % min distance for bonded atoms           
Inputs.randstate     = 0;             
[PP,D,Lmat,Umat]     = setupMC(porig,Inputs);            
[dim,n]              = size(PP);
pars.LOWBD           = full(Lmat);
pars.UPPBD           = full(Umat);
pars.PP              = PP;
D                    = sparse(D);
 

%PREEEDM% ------------------------------------------------------------------
pars.draw  = 1; 
Out        = PREEEDM(D,dim,pars); 






