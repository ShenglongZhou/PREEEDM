function [D, dim, pars] = ArtiMC(s, nf, dim, range, model, anchor)


n       = s^3;
PP      = zeros(dim,n);
[x,y,z] = meshgrid(0:1:s-1, 0:1:s-1, 0:1:s-1);

for i=0:s-1
    for j=0:s-1
    for k=0:s-1            
        PP(:,1+i+s*j+s^2*k)=[x(k+1,j+1,i+1); y(k+1,j+1,i+1); z(k+1,j+1,i+1)];
    end
    end
end

if anchor>0
    PA=(s-1)*rand(3,anchor);
else
    PA=[];
end
m = size(PA,2);I=[];

for i=1:n
    for j=1:m
        if norm(PA(:,j)-PP(:,i))==0; I(j)=i; break; end
    end
end

Ic =  setdiff(1:n,I);
PP = [PA PP(:,Ic)]; 
D  = squareform(pdist(PP'));  
n  = size(D,2); 
 
if model==1
    P = 1+PP(1,:)+s*PP(2,:)+s^2*PP(3,:);
    n = size(PP,2);
    H = zeros(n,n);
    for i=1:n
        for j=i+1:n 
        if  abs( P(i)-P(j) )<=range;
            H(i,j)=1;H(j,i)=1;      
        end
        end
    end
    D           = H.*D;
 
    U           = sqrt(3)*(s-1)*ones(n); %UPPER BOUND
    U(H==1)     = max(D(H==1));
    L           = ones(n);               %LOWER BOUND
    L(1:n+1:n^2)= 0;
else
    D(D>range)  = 0; 
    U           = sqrt(3)*(s-1)*ones(n); %UPPER BOUND
    L           = ones(n);               %LOWER BOUND 
    L(1:n+1:n^2)= 0;  
    pars.range  = range;
end

D         = D.*(abs(1+randn(size(D))*nf));
D         = (D+D')/2; 
D(1:m,1:m)= squareform(pdist(PP(:,1:m)'));
pars.PP   = PP;
pars.m    = m;    
pars.UPPBD= U;
pars.LOWBD= L;

end