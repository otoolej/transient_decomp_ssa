function [hfig]=set_figure(h,CLEAR_FIG)
if(nargin<2 || isempty(CLEAR_FIG)), CLEAR_FIG=1; end


if((ishandle(h) && strcmp(get(h, 'type'), 'figure')))
    set(0,'currentfigure',h);
    if(CLEAR_FIG), clf; end
    hold all;
    hfig=h;
else
    hfig=figure(h); 
    if(CLEAR_FIG), clf; end
    hold all;
end
hfig=gcf;
