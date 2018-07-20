%%*****************************************************************
%% This file is part of DISCO: 
%% Copyright (c) 2009
%% Kim-Chuan Toh
%% Last Modified: 16 Sep 2012
%%*****************************************************************

  function [Dv,Lv,Uv] = noisifydistances(Dv,noiseType,noiseLevel,minLowerBound)
  rng('default');
  rng('shuffle');
  m = length(Dv);

  switch(noiseType)
  case 'normal'
     const = sqrt(pi/2); 
     Lv = max(minLowerBound, Dv.*(1-const*noiseLevel*abs(randn(m,1))));
     Uv = max(minLowerBound, Dv.*(1+const*noiseLevel*abs(randn(m,1))));
  case 'uniform'
     const = 2; 
     Lv = max(minLowerBound, Dv.*(1-const*noiseLevel*rand(m,1)));
     Uv = max(minLowerBound, Dv.*(1+const*noiseLevel*rand(m,1)));
  otherwise
     error('Noise type not implemented yet in noisifydistances.m.');
  end

  Dv = 1/2 * (Lv+Uv);
%%*****************************************************************
