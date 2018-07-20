%%*****************************************************************
%% This file is part of DISCO: 
%% Copyright (c) 2009
%% Kim-Chuan Toh
%% Last Modified: 16 Sep 2012
%%*****************************************************************

  function [Di,Dj,Dv] = prunedistances(Di,Dj,Dv,nAtoms,sparseLevel)

  nEdgesBefore = length(Di);

  [T,D] = computespanningtree(sparse(Di,Dj,Dv,nAtoms,nAtoms));

  [Ti,Tj,Tv] = find(T);
  [Di,Dj,Dv] = find(D);
  sparseLevelNew = (sparseLevel * nEdgesBefore - nAtoms + 1) / ...
                   (nEdgesBefore - nAtoms + 1);
  Dv = Dv .* (rand(length(Dv),1) < sparseLevelNew);
%%
  Di = [Di; Ti];
  Dj = [Dj; Tj];
  Dv = [Dv; Tv];
%%
  idx = find(Dv > 0); 
  Dv = Dv(idx); 
  Di = Di(idx); 
  Dj = Dj(idx); 
  nEdgesAfter = length(Dv); 
  fprintf('prunedistances: Proportion of edges given within the cut-off radius = %3.2f%%\n',...
           nEdgesAfter/nEdgesBefore*100);
  fprintf('max, min distance = %3.2f, %3.2f\n',max(Dv),min(Dv)); 
%%****************************************************************
