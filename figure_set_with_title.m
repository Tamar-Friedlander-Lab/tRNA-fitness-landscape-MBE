%Set default Figure/Axe/object Properties to my liking
%Axes color
cColor = [0.5 0.5 0.5];%dark grey
%Set default Figure properties
%set(0, 'DefaultPaperPositionMode', 'auto'); %Use this if you want the print to be
%same format as screen
%Change these properties for setting the page layout
margin = 4;
width = 25;
set(0,'DefaultFigurePaperUnits','centimeters')
set(0,'DefaultFigurePaperPosition',[margin margin width 0.61*width])%a figure fits in
%a page, approx golden ratio
%Set default Axe properties
% set(0,'DefaultAxesColor',cColor)
%set(0,'DefaultAxesFontName','Helvetica')
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultAxesFontsize',12)
set(0,'DefaultAxesLineWidth',1);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesTickDir','out')
% set(0,'DefaultAxesXcolor',cColor)
% set(0,'DefaultAxesYcolor',cColor)
% set(0,'DefaultAxesZcolor',cColor)
set(0,'DefaultAxesBox','off')
set(0,'DefaultAxesTickLength',[.01 .01])
set(0,'DefaultAxesXMinorTick','off')
set(0,'DefaultAxesYMinorTick','off')
%left = 0.12; bottom = 0.12; width = 0.78; height = 0.78;
% left = 0.19; bottom = 0.19; width = 0.78; height = 0.78; 
left = 0.19; bottom = 0.12; width = 0.78; height = 0.78; % 0.78
set(0,'DefaultAxesPosition',[left bottom width height]); % 16/3/2015: changed the original values. Otherwise the xlabel is trimmed. 
%set(0,'DefaultAxesPosition',[0.1 0.1 0.8 0.8])
%Set default Text properties
%set(0,'DefaultTextColor',cColor)
%set(0,'DefaultTextFontName','Helvetica')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultTextFontSize',12)
set(0,'DefaultTextFontWeight','bold')
set(0,'DefaultTextInterpreter','Tex') 
%set(0,'DefaultTextInterpreter','Latex') % Use only if you want special characters in
%texts. You must find a way to change the font to helvetica though or worse idea,
%change all other fonts
%Set default Line properties
set(0,'DefaultLineLineWidth',2)
clear cColor;
clear margin;
clear width;