%%**********************************************************************
%% This file is part of DISCO: 
%% Copyright (c) 2009
%% Kim-Chuan Toh
%% Last Modified: 16 Sep 2012
%%**********************************************************************
%%**********************************************************************
%% The inputs are:
%%   o noiseLevel
%%   o noiseType
%%   o radius
%%   o sparseLevel
%%   o minLowerBound
%%
%% The outputs are:
%%   o A
%%   o Dmat, Lmat, Umat
%%
%% Please note:
%% This function reads the PDB file given in [filename]
%% [A,Dmat,Lmat,Umat,nLooseAtoms] = setup(filename,Inputs)
%%**********************************************************************

   function [PP,Dmat,Lmat,Umat,nLooseAtoms] = setupMC(PP,Inputs)

%% Read inputs.

   randstate = 0; 
   if isfield(Inputs,'randstate'); randstate = Inputs.randstate; end
   rand('seed',randstate);
   randn('seed',randstate);  

   noiseLevel  = Inputs.noiseLevel;
   noiseType   = Inputs.noiseType;
   sparseLevel = Inputs.sparseLevel;
   radius      = Inputs.radius;
   minLowerBound = Inputs.minLowerBound;

%% Read A, center A at origin, permute the columns of A.

   timeTaken = cputime;
   nDim   = size(PP,1); 
   nAtoms = size(PP,2);    
   centroidA = PP*ones(nAtoms,1)/nAtoms; 
   PP = PP - centroidA*ones(1,nAtoms);
   PP = PP(:,randperm(nAtoms));  

%% Compute D, L, U.
   [Di,Dj,Dv] = computedistances(PP,radius);
   [Di,Dj,Dv] = prunedistances(Di,Dj,Dv,nAtoms,sparseLevel);
   [Dv,Lv,Uv] = noisifydistances(Dv,noiseType,noiseLevel,minLowerBound);

%% Setup oD, oL, oU.
   Dmat = sparse(Di,Dj,Dv,nAtoms,nAtoms); 
   Lmat = sparse(Di,Dj,Lv,nAtoms,nAtoms); 
   Umat = sparse(Di,Dj,Uv,nAtoms,nAtoms); 
   degree = sum(spones(Dmat + Dmat'));
   nLooseAtoms = length(find(degree < nDim+1)); 
  
   Dmat = triu(Dmat,1)+triu(Dmat,1)'; 
   Lmat = triu(Lmat,1)+triu(Lmat,1)'; 
   Umat = triu(Umat,1)+triu(Umat,1)';
%    H    = spones(Dmat);
%    Umat = Umat+(1-H)*inf;

%% Print
   timeTaken = cputime - timeTaken;
 
   fprintf('Time taken in setup = %3.1f s\n',timeTaken);
   fprintf('Number of atoms with < 4 neighbors = %2.0d(%3.2f%%)\n',...
            nLooseAtoms,nLooseAtoms/nAtoms*100); 
%%****************************************************************
