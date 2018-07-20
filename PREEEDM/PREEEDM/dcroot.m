
function x =dcroot( w, a, Ind, d2, lbd, ubd  )

%      min    0.5*(x-w)^2 + a*|sqrt(x)-sqrt(d2)|  
%      s.t.   lbd <= x <= ubd.
%
% where a>=0, Ind=find(a>0)

    x  = min( ubd, max( lbd, w ) ); 
    I1 = Ind(find(d2(Ind)<=lbd(Ind)));  
    if ~ isempty(I1)        
         x(I1) = dcroot_plus(w(I1),a(I1),lbd(I1),ubd(I1));
    end
    I2 = Ind(find(d2(Ind)>=ubd(Ind))); 
    if ~ isempty(I2)        
         x(I2) = dcroot_minus(w(I2),a(I2),lbd(I2),ubd(I2));
    end
    I3 = setdiff(Ind, union(I1,I2)); 
    if ~ isempty(I3)        
         x(I3) = dcroot_plus_minus(w(I3),a(I3),d2(I3),lbd(I3),ubd(I3));
    end
    x  = real(x);
end


function x = dcroot_plus( w, a, lbd, ubd )
% min 0.5*(x-w)^2+a*sqrt(x) s.t. lbd<=x<=ubd. (where a>0)
    x  = ubd;
    a2 = a/2;
    l1 = lbd + a2./sqrt(lbd);
    u1 = ubd + a2./sqrt(ubd);
    I1 = find(w<=l1);
    if ~ isempty(I1)        
         x(I1)=lbd(I1);
    end
    I2 = find(w>l1 & w<u1);
    if ~ isempty(I2)   
         v   = w(I2)/3;
         x(I2)= 2*v.*(1+cos((2*pi/3)-(2/3)*acos(a2(I2)./(v.^1.5)/2)));
    end
end


function x = dcroot_minus( w, a, lbd, ubd )
% min 0.5*(x-w)^2-a*sqrt(x) s.t. lbd<=x<=ubd. (where a>0)
    x  = lbd;
    a2 = a/2;
    l1 = lbd - a2./sqrt(lbd);
    u1 = ubd - a2./sqrt(ubd);
    I1 = find(w>=u1);
    if ~ isempty(I1)        
         x(I1)=ubd(I1);
    end
    I2 = find(w>l1 & w<u1);
    if ~ isempty(I2)   
         x(I2)= dcroot_minus_med(w(I2),a(I2));
    end
end

function x = dcroot_minus_med( w, a)
% min 0.5*(x-w)^2-a*sqrt(x) s.t. x>=0. (where a>0)
    x  = w;
    a  = a/4;
    w  = w/3;
    a2 = a.^2; 
    w3 = w.^3; 
    d  = a2-w3;
    I1 = find(d<0); 
    if ~ isempty(I1) 
         x(I1) = (2*w(I1)).*(1+cos((2/3)*acos(sqrt(a2(I1)./w3(I1)))));
    end    
    I2 = find( d>=0 & w>=0);
    if ~ isempty(I2)
         x(I2) = (a(I2)+sqrt(d(I2))).^(1/3)+(a(I2)-sqrt(d(I2))).^(1/3);
         x(I2) = x(I2).*x(I2);
    end
    I3 = find( d>=0 & w<0);
    if ~ isempty(I3)
         x(I3) = (a(I3)+sqrt(d(I3))).^(1/3)-(sqrt(d(I3))-a(I3)).^(1/3);
         x(I3) = x(I3).*x(I3);
    end
end


function x = dcroot_plus_minus( w, a, d2, lbd, ubd )
% min 0.5*(x-w)^2+a*|sqrt(x)-sqrt(d2)| s.t. lbd<=x<=ubd. 
% where a>0, lbd<d2<ubd)
    x  = ubd;
    a2 = a/2;
    l1 = lbd - a2./sqrt(lbd);
    dl = d2  - a2./sqrt(d2);
    du = d2  + a2./sqrt(d2);
    u1 = ubd + a2./sqrt(ubd);
    I1 = find(w<=l1);
    if ~ isempty(I1)        
         x(I1)=lbd(I1);
    end
    I2 = find(w>l1 & w<dl);
    if ~ isempty(I2)   
         x(I2) = dcroot_minus_med(w(I2), a(I2));
    end
    I3 = find(w>=dl & w<=du);
    if ~ isempty(I3)   
         x(I3) = d2(I3);
    end
    I4 = find(w>du & w<u1);
    if ~ isempty(I4)   
         v     = w(I4)/3;      
         x(I4) = 2*v.*(1+cos((2*pi/3)-(2/3)*acos(a2(I4)./(v.^1.5)/2)));
    end
end

