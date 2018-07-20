function [D, dim, pars] = EDMExamples( n0,m, nf, range, draw)  
% This file was modified from the code of EDME which can be download at
% http://www.personal.soton.ac.uk/hdqi/papers.html
% It generates the letter E D M in the graph
% Inputs:
%    n0    -- number of sensors in total, known and unknown
%    m     -- number of known sensors (anchors)
%    range -- communicating radio range
%    draw  -- =1 draw the groundtruth graph, =0 otherwise
while 1
	seed       = 1;
    [X,idxE,idxD,idxM] = genEDM(n0,seed,draw);
    n          = size(X,1);
    D          = squareform(pdist(X));
    Dold       = D;
    D(D>range) = 0;
    rng('default');
    rng('shuffle');
    D          = triu(D).*abs(1+nf*randn(n));
    D          = D+D';
    D(1:m,1:m) = Dold(1:m,1:m);
    
    [~,bins]   = graphconncomp(sparse(D));
    if all(bins == 1);   break;
    else  error('Graph is not connected, try to increase range...')
    end
end
pars.m     = m;
pars.n     = n;
pars.PP    = X';     
dim        = 2;
pars.range = range;
pars.E     = idxE;
pars.D     = idxD;
pars.M     = idxM;
D          = sparse(D);
end  

function [X,idxE,idxD,idxM] = genEDM(n0,seed,draw)
if nargin < 3
    draw = 0;
end
n=ceil(n0/0.7);
rng('default');
rng(seed);
n = 3*ceil(n/3);
x = 115*rand(n,1);
y = 50*rand(n,1);

for i = 1:n
    if x(i)>10 && x(i)<40 && ((y(i)>10 && y(i)<20) || (y(i)>30 && y(i)<40))
        x(i) = 0;
        y(i) = 0;
    end
end

for i = 1:n
    if x(i)>50 && x(i)<75 && (norm([x(i),y(i)]-[50,25])<15 || norm([x(i),y(i)]-[50,25])>25)
        x(i) = 0;
        y(i) = 0;
    end
end

for i = 1:n
    if x(i)>85 && x(i)<95 && (y(i)<2.5*(95-x(i)) || y(i)>2.5*(105-x(i)))
        x(i) = 0;
        y(i) = 0;
    end
end

for i = 1:n
    if x(i)>95 && x(i)<105 && (y(i)<2.5*(x(i)-95) || y(i)>2.5*(x(i)-85))
        x(i) = 0;
        y(i) = 0;
    end
end

X = [x,y];
k = 1;
for i = 1:n
    if X(i,1)~=0 && X(i,2)~=0
        Y(k,:) = X(i,:);
        k = k+1;
    end
end
 
k =k-1;  clear X;
if k>n0
    I = randperm(k);
    X = Y(I(1:n0),:);
%     X = Y(k-n0+1:k,:); 
elseif k<n0
    I  = randperm(k); 
    X0 = Y(I(1:n0-k),:)+rand(size(Y(I(1:n0-k),:)));
    X  = [Y; X0];
else
    X = Y;
end

X(1,:)=[110,8];
mX    = max(X(:));
X     = X/mX;
xp    = X(:,1);
idxE  = find(xp<40/mX);
idxD  = find(xp>=40/mX & xp<75/mX);
idxM  = find(xp>=75/mX);
    
    
if draw
    figure;
    scatter(X(idxE,1),X(idxE,2),12,'filled');
    hold on;
    scatter(X(idxD,1),X(idxD,2),12,'filled');
    scatter(X(idxM,1),X(idxM,2),12,'filled');
    set(gca,'FontName','Times','FontSize',8);
    xlabel(['Ground truth EDMnetwork']);
    axis equal;
    axis([-5/mX 120/mX -5/mX 55/mX]);
end

end
