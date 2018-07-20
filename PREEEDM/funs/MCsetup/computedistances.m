%%*****************************************************************
%% This file is part of DISCO: 
%% Copyright (c) 2009
%% Kim-Chuan Toh
%% Last Modified: 16 Sep 2012
%%*****************************************************************
%%*****************************************************************
%% find all the pairs of atoms with distances < radius. 
%%
%%*****************************************************************

  function [Di,Dj,Dv] = computedistances(A,radius)

  [nDimensions,nAtoms] = size(A);

  Di = [];
  Dj = [];
  Dv = [];
  for i = 1:nAtoms
     m = nAtoms - i;
     AiAj = A * ...
           (sparse(i,1:m,1,nAtoms,m) - sparse((i+1):nAtoms,1:m,1,nAtoms,m));
     d = sqrt(sum(AiAj.^2))';
     idx = find(d < radius); 
     Di = [Di; i * ones(length(idx),1)];
     Dj = [Dj; i + idx];
     Dv = [Dv; d(idx)];
  end
%%*****************************************************************
