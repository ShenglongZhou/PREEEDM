clc; clear; close all;    

%Initialization% ----------------------------------------------------------
s         = 8;
nf        = 0.1; 
model     = 2;
anchor    = 0;
if model==1; R = s^2;
else         R = 2;
end
[D, dim, pars] = ArtiMC(s, nf, 3, R, model, anchor);   

%PREEEDM% ------------------------------------------------------------------

pars.draw  = 1; 
Out        = PREEEDM(D,dim,pars) ; 
 
