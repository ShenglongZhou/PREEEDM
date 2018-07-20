function  plot_SNL(X,RMSD,refine,pars,CM)
% X      -- dim-by-(n-m) matrix containing coordinates of (n-m) points, where m>=0
% RMSD   -- Root Mean Square Distance (RMSD) before refinement if refine=0
%           Root Mean Square Distance (RMSD) after refinement  if refine=1
% refine -- 0 plot X before refinement
%           1 plot X before refinement
% pars   -- pars.PP:  dim-by-n matrix of coordinates of n points with first
%                     m(>=0) columns being anchors
%           pars.m:   anchors' number
%           pars.E:   coordinates of points forming shap 'E' in EDM network
%           pars.D:   coordinates of points forming shap 'D' in EDM network
%           pars.M:   coordinates of points forming shap 'M' in EDM network
% CM     -- Connectivity matrix of a graph, i.e.,
%                     Dis = squareform(pdist(pars.PP'));
%                     mD  = max(Dis(:));
%                     CM  = full( spones( Dis(Dis<0.1*mD) ) );

if  ~isfield(pars,'PP'); return; end 

if ~isfield(pars,'E')  
    
    if isfield(pars, 'm'); m=pars.m; else m=0;  end 
    PP   = pars.PP;    
    [r,n]= size(PP);
    Ta   = 1:m;
    Ts   = m+1:n;
    if r==2  % plot in 2 dimensional space
            if nargin==5
               for i=1:n
               for j=1:n
               if  CM(i,j)==1;  
                   line([PP(1,i) PP(1,j)], [PP(2,i) PP(2,j)],...
                        'linewidth',.1,'color','b'); hold on
               end
               end
               end
            end
            Af=[PP(:,Ta) X]; 
            set(gca,'FontName','Times','FontSize',8)           
            plot(PP(1,Ts),PP(2,Ts),'bo','markersize',4.5);hold on
            plot(X(1,:),X(2,:),'m*','markersize',3);hold on   
            plot(PP(1,Ta),PP(2,Ta),'gs','markersize',4,'linewidth',2);     
            for i=1:n-m; 
                line([X(1,i) PP(1,m+i)], [X(2,i) PP(2,m+i)],...
                    'linewidth',.1,'color','b'); hold on
            end                   
            plot(PP(1,Ta),PP(2,Ta),'gs','markersize',4,'linewidth',2);
            ZZ = [PP Af]';
            if ~isempty(RMSD)
            if refine==1; title(['rRMSD = ', sprintf('%4.2e', RMSD)],'FontName','Times','FontSize',8);  
            else          title(['RMSD = ',  sprintf('%4.2e', RMSD)],'FontName','Times','FontSize',8); 
            end
            end
            axis([min(min(ZZ(:,1))) max(max(ZZ(:,1))) min(min(ZZ(:,2))) max(max(ZZ(:,2)))]) 
            hold on
    else    % plot in 3 dimensional space
        if nargin==5
            set(gca,'FontName','Times','FontSize',8)           
            plot3(PP(1,Ts),PP(2,Ts),PP(3,Ts),'bo','markersize',4.5);hold on          
               for i=1:n
               for j=1:n
               if CM(i,j)==1;  
                   line([PP(1,i) PP(1,j)], [PP(2,i) PP(2,j)], [PP(3,i) PP(3,j)],...
                        'linewidth',.1,'color','b'); hold on
               end
               end
               end
        end
            
        Af=[PP(:,Ta) X];   
            set(gca,'FontName','Times','FontSize',8) 
            plot3(X(1,:),X(2,:),X(3,:),'m*','markersize',3);hold on 
            for i=1:n-m; 
                line([X(1,i) PP(1,m+i)], [X(2,i) PP(2,m+i)], [X(3,i) PP(3,m+i)],...
                    'linewidth',.1,'color','m'); hold on
            end
            plot3(PP(1,Ts),PP(2,Ts),PP(3,Ts),'bo','markersize',4.5);hold on
            plot3(X(1,:),X(2,:),X(3,:),'m*','markersize',3);hold on   
%            plot3(PP(1,T),PP(2,T),PP(3,T),'gs','markersize',4,'linewidth',2); hold on            
            ZZ = [PP Af]';
            plot3(PP(1,Ta),PP(2,Ta),PP(3,Ta),'gs','markersize',4,'linewidth',2); hold on  
            if ~isempty(RMSD)
            if refine==1; title(['rRMSD = ', sprintf('%4.2e', RMSD)],'FontName','Times','FontSize',8);  
            else          title(['RMSD = ',  sprintf('%4.2e', RMSD)],'FontName','Times','FontSize',8); 
            end
            end
            grid on;
            axis([min(min(ZZ(:,1))) max(max(ZZ(:,1))) min(min(ZZ(:,2))) max(max(ZZ(:,2)))...
                min(min(ZZ(:,3))) max(max(ZZ(:,3)))]) 
            hold on             
    end    
else       % plot EDM in 2 dimensional space
    PX=[pars.PP(:,1:pars.m) X];
    scatter(PX(1,pars.E),PX(2,pars.E),12,'filled');
	hold on
	scatter(PX(1,pars.D),PX(2,pars.D),12,'filled');
	scatter(PX(1,pars.M),PX(2,pars.M),12,'filled');
 	scatter(PX(1,1:pars.m),PX(2,1:pars.m),12,'filled','k');
    set(gca,'FontName','Times','FontSize',8);
    if refine==1; title(['rRMSD = ', sprintf('%4.2e', RMSD)],'FontName','Times','FontSize',8);  
    else          title(['RMSD = ',  sprintf('%4.2e', RMSD)],'FontName','Times','FontSize',8); 
    end
	axis equal;
	axis([-0.05*max(PX(1,:)) max(PX(1,:))*1.05...  
          -0.1*max(PX(2,:)) max(PX(2,:))*1.075 ]);
end
    
end




