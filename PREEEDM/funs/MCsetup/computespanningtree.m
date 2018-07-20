%%*****************************************************************
%% This file is part of DISCO: 
%% Copyright (c) 2009
%% Kim-Chuan Toh
%% Last Modified: 16 Sep 2012
%%*****************************************************************
%% computespanningtree
%%
%% Assumption: D is an upper-triangular distance matrix.
%%
%% Running time: O(|E|), where E is the edge set.
%%
%% [T,D] = computespanningtree(D)
%%*****************************************************************

function [T,D] = computespanningtree(D)

% i(k) and j(k) are connected by the k-th edge.
[i,j]  = find(D);

% Permute edges.
nEdges = length(i);
p      = randperm(nEdges);
i      = i(p);
j      = j(p);

nVertices = size(D,1);
T         = sparse(nVertices,nVertices);
group     = zeros(nVertices,1);
nGroups   = 0;
for iEdge = 1:nEdges
    iGroup = group(i(iEdge));
    jGroup = group(j(iEdge));
    if iGroup == 0 && jGroup == 0
        T(i(iEdge),j(iEdge)) = 1;
        nGroups              = nGroups + 1;
        group(i(iEdge))      = nGroups;
        group(j(iEdge))      = nGroups;
    elseif iGroup == 0 && jGroup ~= 0
        T(i(iEdge),j(iEdge)) = 1;
        group(i(iEdge))      = jGroup;
    elseif iGroup ~= 0 && jGroup == 0
        T(i(iEdge),j(iEdge)) = 1;
        group(j(iEdge))      = group(i(iEdge));
    else % iGroup ~= 0 && jGroup ~= 0
        if iGroup == jGroup
            % Do nothing.
        else
            T(i(iEdge),j(iEdge))   = 1;
            group(group == jGroup) = iGroup;
        end
    end
end

T = T .* D;
for iVertex = 1:nVertices
    D(:,iVertex) = (~T(:,iVertex)) .* D(:,iVertex);
end

% The following code doesn't work cos ~T can be too big:
%     D = (~T) .* D;
% So do it column by column.

