function [] = plot_ssa_components(x_components, U, fig_num)
if(nargin < 3 || isempty(fig_num)), fig_num = 1; end


[L, N] = size(x_components);
x = sum(x_components);


set_figure(fig_num);

%---------------------------------------------------------------------
% 1. plot the components
%---------------------------------------------------------------------
hax = subplot(4,3,[1,2,4,5,7,8]); 
hold all;
plot_components_(x_components);
hax.YAxis.Visible = 'off';
title('components');

subplot(4,3,[10,11]); hold all;
plot(x); xlim([1 N]);
xlabel('samples');

%---------------------------------------------------------------------
% 2. plot frequency response of filters:
%---------------------------------------------------------------------
L_pad=1024; 
if(L>L_pad)
    L_pad=L; 
end
L_pad_h=ceil(L_pad/2);

for n=1:L
    H=fft(U(:,n),L_pad);
    His(n,:)=abs(H(1:L_pad_h)).^2;
end

hax = subplot(4,3,[3,6,9]); 
hold all;
plot_components_(His,[],linspace(0,0.5,L_pad_h));
hax.YAxis.Visible = 'off';
grid('on');
% xlabel('normalised frequency');
title('eigenvectors');

%---------------------------------------------------------------------
% 3. plot the log-log spectrum
%---------------------------------------------------------------------
subplot(4,3,12); 
hold all;
X=fft(x,L_pad); pxx=abs(X(1:L_pad_h)).^2;
pspec(x);




function [hlines]=plot_components_(x_comp,y_gap,xticks)
%---------------------------------------------------------------------
% plot all components in one axis
%---------------------------------------------------------------------
if(nargin<2 || isempty(y_gap)), y_gap=0; end
if(nargin<3 || isempty(xticks)), xticks=[]; end


[L,N]=size(x_comp);
if(isempty(xticks))
    xticks=1:N;
end


y_gap=mean(std(x_comp'));

hlines=zeros(1,L);
all_gap=0;
for l=1:L
    if(l>1)
        yheight=max(x_comp(l,:))-mean(x_comp(l,:));
        all_gap=y_gap+abs(yl)+yheight;
    else
        x_comp(1,:)=x_comp(1,:)-nanmean(x_comp(1,:));
    end
    
    hlines(l)=plot(xticks,x_comp(l,:)-all_gap);
    yl=min(get(hlines(l),'ydata'));
end


xlim([xticks(1) xticks(end)]);
ylim([min(get(hlines(end),'ydata')) max(get(hlines(1),'ydata'))]);



function h_freq=pspec(x,fs,fignum)
%---------------------------------------------------------------------
% log-log plot of spectrum
%---------------------------------------------------------------------
if(nargin<2) fs=1; end
if(nargin<3) fignum=[]; end


N=length(x); Nh=ceil(N/2);
if(~isempty(fignum))
    figure(fignum); clf;
end
c=fft(x); ch=c(1:Nh); cpow=abs(ch).^2;


%---------------------------------------------------------------------
% Calculating dB values?
%---------------------------------------------------------------------
cpow=20*log10( cpow+eps );

k=0:Nh-1; k=k+1;
f=linspace(0,fs/2,Nh);

%---------------------------------------------------------------------
% Do the plot
%---------------------------------------------------------------------
h_freq=plot(f,cpow(k)); 


axis('tight'); grid('on');
ylabel('magnitude (dB)');
xlabel('normalised frequency');