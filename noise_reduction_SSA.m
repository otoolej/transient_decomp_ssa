function [yb,d]=noise_reduction_SSA(x,L,method,use_dct,dct_keep,remove_mean,d_max,db_plot)
if(nargin<2 || isempty(L)), L=20; end
if(nargin<3 || isempty(method)), method='Vautard_Ghil'; end
if(nargin<4 || isempty(use_dct)), use_dct=false; end
if(nargin<5 || isempty(dct_keep)), dct_keep=1; end
if(nargin<6 || isempty(remove_mean)), remove_mean=false; end
if(nargin<7 || isempty(d_max)), d_max=[]; end
if(nargin<8 || isempty(db_plot)), db_plot=false; end

DBverbose=0;

N=length(x);
x=x(:)';

SURR_WGN=1;

%---------------------------------------------------------------------
% remove DC?
%---------------------------------------------------------------------
if(remove_mean)
    xmean=nanmean(x);
    x=x-xmean;
end


%---------------------------------------------------------------------
% operate on DCT, not the time-domain signal; equavalent to 90 degree
% rotation in the TF domain
%---------------------------------------------------------------------
if(use_dct)
    x_time=x;
    x=dct(x);
    x_dct=x;
    if(dct_keep~=1)
        L_cut=round(length(x).*dct_keep);
        x((L_cut+1):end)=0;
    end
end

%---------------------------------------------------------------------
% set which SSA to use (either diagonal averaging or forward-reverse 
% filtering)
%---------------------------------------------------------------------
ssa=@(x,L,weights) ssa_filter_bank_approach(x,L,weights);





switch lower(method)
  case 'vautard_ghil'
    %---------------------------------------------------------------------
    % null hypothesis testing using surrogate data (from [4], above)
    %---------------------------------------------------------------------
    N_iter=100;
    
    
    if(SURR_WGN)
        w=randn(N_iter,N);
    else
        w=zeros(N_iter,N);
        for n=1:N_iter
            w(n,:)=dct(x_time(randperm(N)));
        end
    end

    d=0; p=0;
    
    if(isempty(d_max))
        d_max=L;
    end
    
    while(d==0 && p<=d_max)
        weights=ones(1,L); weights(1:p)=0;
        y=do_ssa(x,L,weights,ssa);
        
        
        % generate auto-correlation function on signal:
        [ry,l]=xcorr(y,L);
        izerolag=find(l==0);
        ry=ry(izerolag:(end-1));

        % and same for surrogate data:
        rym=NaN(N_iter,L);
        for n=1:N_iter
            rym(n,:)=ssa_xcorr(w(n,:),L,weights,izerolag,ssa);
        end            
        % ci_upper=nanmean(rym) + 1.96.*nanstd(rym);
        % ci_lower=nanmean(rym) - 1.96.*nanstd(rym);
        
        % 99 percentile intervals:
        ci_upper=nanmean(rym) + 2.576.*nanstd(rym);
        ci_lower=nanmean(rym) - 2.576.*nanstd(rym);
        
        
        iry_outofbounds_upper=find(ry>ci_upper);
        iry_outofbounds_lower=find(ry<ci_lower);
           
        b1=NaN; b2=NaN;
        if(~isempty(iry_outofbounds_upper))
            b1=max([ry(iry_outofbounds_upper)./ci_upper(iry_outofbounds_upper)]);
        end
        if(~isempty(iry_outofbounds_lower))
            b2=max([ry(iry_outofbounds_lower)./ci_lower(iry_outofbounds_lower)]);
        end
        bbeta=sqrt(max([b1 b2]));
        if(isnan(bbeta))
            bbeta=1;
        end
        
        scaled_ci_upper=(bbeta.^2).*ci_upper;
        scaled_ci_lower=(bbeta.^2).*ci_lower;
        iry_outofbounds_upper=find(ry>scaled_ci_upper);
        iry_outofbounds_lower=find(ry<scaled_ci_lower);
        
        if((isempty(iry_outofbounds_lower) && isempty(iry_outofbounds_upper)))
            d=p;
        end
        
        
        DBplot_iter=0;
        if(DBplot_iter)
            set_figure(4); 
            dispVars(p);
            plot(ci_upper,'k--+'); plot(ci_lower,'k--+'); plot(nanmean(rym),'b-+');
            plot(ry,'ro','linewidth',2);
            plot(scaled_ci_upper,'g--+'); plot(scaled_ci_lower,'c--+');
            disp('--- paused; hit key to continue ---'); pause;
        end
        

        p=p+1;
    end
    if(d==0), d=p-1; end
    
    if(DBverbose)
        fprintf('¬¬¬ VG END; d=%d\n',d);
    end
    


    % remove noise components:
    weights=zeros(1,L); weights(1:d)=1;
    yb=do_ssa(x,L,weights,ssa);

    
  case 'celka_mdl'    
    %---------------------------------------------------------------------
    % Minimum descriptive length ([3] and [4] above)
    %---------------------------------------------------------------------
    [y,~,eigv]=ssa(x,L,[]);
    
    
    K=N-L+1;    
    mdl=NaN(1,L-1);
    for n=1:L-1
        mdl(n)=min_description_length(eigv,n,K);
    end
    
    [~,imin]=min(mdl);
    
    if(~isempty(d_max))
        if(imin>d_max)
            imin=d_max;
        end
    end
    

    DBmdl=0;
    if(DBmdl)
        set_figure(3);
        plot(mdl,'-o');
        stem(imin,mdl(imin),'r+');
    end
    

    d=imin;
    weights=zeros(1,L); weights(1:d)=1;
    yb=do_ssa(x,L,weights,ssa);
    if(DBverbose)
        fprintf('¬¬¬ MDL END; d=%d\n',d);
    end
    

    
    
  case 'eigenvalue_ranking'    
    %---------------------------------------------------------------------
    % typical approach:
    % rank eigenvalues and select the first X by finding largest jump 
    % ranked eigenvalues (using the derivative)
    %---------------------------------------------------------------------
    [~,~,eigv]=ssa(x,L,[]);    
    [~,ip]=max(abs(diff(log(eigv))));
    
    if(~isempty(d_max))
        if(ip>d_max)
            ip=d_max;
        end
    end

    
    weights=zeros(1,L); weights(1:ip)=1;
    yb=do_ssa(x,L,weights,ssa);

    
    DBplot_eigenvalues=0;
    if(DBplot_eigenvalues)
        set_figure(12);
        subplot(4,1,[1:3]);
        plot(log(eigv),'-o');
        subplot(4,1,4);
        plot(abs(diff(log(eigv))),'-o');
        fprintf('Eigenvalue ranking method: use index=%d\n',ip);        
    end
    if(DBverbose)
        fprintf('¬¬¬ EVR end; d=%d\n',ip);
    end
    d=ip;
    
  otherwise
    error('method??');
end




if(use_dct)
    yb=idct(yb);
end


if(remove_mean)
    yb=yb+xmean;
end
   

if(db_plot)
    set_figure(2);
    if(use_dct), x=idct(x_dct); end
    
    hx(1)=subplot(3,1,1); hold all;
    plot(x);
    title('signal');
    hx(2)=subplot(3,1,2); hold all;
    plot(yb);
    title('component');
    hx(3)=subplot(3,1,3); hold all;
    plot(x-yb); 
    title('signal minus component');
    linkaxes(hx,'x');
    axis('tight');
    xlabel('samples');
    drawnow;
end



function mdl=min_description_length(eigv,k,Nt)
%---------------------------------------------------------------------
% estimate Rissanen's MDL
%---------------------------------------------------------------------
L=length(eigv);

ik_last=(k+1):L;
ik_first=1:k;

nf=L*k-(k^2/2)+(k/2)+1;

A=prod( eigv(ik_last).^(1/(L-k)) );
B=(1/(L-k))*sum( eigv(ik_last) );
C=sum( log(eigv(ik_first).*sqrt(2/Nt)) );

mdl=-log((A/B)^(Nt*(L-k))) + nf*(0.5+log(32)) - (nf/k)*C;




function y=do_ssa(x,L,weights,ssa)
%---------------------------------------------------------------------
% wrapper for SSA, ignore if all weights =1;
%---------------------------------------------------------------------
if(all(weights==1))
    y=x;
else
    y=ssa(x,L,weights);
end




function rym=ssa_xcorr(w,L,weights,izerolag,ssa)
%---------------------------------------------------------------------
% xcorr on the SSA decomposition
%---------------------------------------------------------------------
yw=ssa(w,L,weights);

rymt=xcorr(yw,L);
rym=rymt(izerolag:(end-1));
