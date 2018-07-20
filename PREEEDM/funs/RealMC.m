function [D, dim, pars] = RealMC(PP, range, sparseLevel, nf )

[dim,n]    = size(PP);

D0         = squareform(pdist(PP'));
tD0        = tril(D0,0);
T          = find(tD0<range & tD0>0);

while 1
I0         = randperm(nnz(T)); 
I          = T(I0(1:floor(sparseLevel*nnz(T)/2)));
LowMat     = zeros(n);
UppMat     = zeros(n);
LowMat(I)  = max(1, (1- abs(nf*randn(size(I)))).*D0(I));
UppMat(I)  = (1+abs(nf*randn(size(I)))).*D0(I);

LowMat     = LowMat + LowMat';
LowMat(1:n+1:n^2)=0;
UppMat     = UppMat + UppMat';
D          = (LowMat+ UppMat)/2;
UppMat(D==0) = inf;
 
pars.LowMat= LowMat;
pars.UppMat= UppMat;
pars.PP    = PP;
D          = sparse(D);
[~,bins]   = graphconncomp(D);
if all(bins== 1); break; else disp('Graph is not collected...'); end

end

end