clc; clear; close all;    

%Initialization% ----------------------------------------------------------
pro      = 1;
n        = 500;  
nf       = 0.1;  

%Examples% ----------------------------------------------------------------
problem   = {'BLTWY06_Inner', 'Tseng07','EDM'} ; 
noisetype = 'multiplicative';

if pro<3
    R = 0.2;   
    m = 5; 
    [D, dim, pars] = SNLExamples(problem{pro}, noisetype, n,m, nf, R);
else   
    R = 0.1;   
    m = 10;    
    [D, dim, pars] = EDMExamples( n, m, nf, R,1);
end
 

%PREEEDM% -----------------------------------------------------------------
pars.draw = 1; 
Out       = PREEEDM(D,dim,pars) ;

