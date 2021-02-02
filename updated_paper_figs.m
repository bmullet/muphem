%%%%%%%%%%%%%%%%%%%%%
%%% Paper Figures %%%
%%%%%%%%%%%%%%%%%%%%%


%% Set default figure paarmeters
% Source: https://dgleich.wordpress.com/2013/06/04/creating-high-quality-graphics-in-matlab-for-papers-and-presentations/

% Defaults for this blog post
width = 4;     % Width in inches
height = 4;    % Height in inches
alw = 1.5;    % AxesLineWidth
fsz = 15;      % Fontsize
lw = 4;      % LineWidth
msz = 8;       % MarkerSize

% The new defaults will not take effect if there are any open figures. To
% use them, we close all figures, and then repeat the first example.
close all;

% The properties we've been using in the figures
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz

% Set the default Size for display
defpos = get(0,'defaultFigurePosition');
set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*100]);

% Set the defaults for saving/printing to a file
set(0,'defaultFigureInvertHardcopy','on'); % This is the default anyway
set(0,'defaultFigurePaperUnits','inches'); % This is the default anyway
defsize = get(gcf, 'PaperSize');
left = (defsize(1)- width)/2;
bottom = (defsize(2)- height)/2;
defsize = [left, bottom, width, height];
set(0, 'defaultFigurePaperPosition', defsize);
set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultTextFontname', 'Arial')

%% Stress regimes
figure
subplot(121)
clrs = parula(3);
xlim([0,2])
ylim([0,2]);
xl = xlim;
plot(xl,[1,1],'Color', clrs(1,:)); hold on;
plot(xl,xl,'Color', clrs(2,:))
plot(xl,2*xl-1,'Color', clrs(3,:))
grid on;
xlabel("k")
ylabel("p/\sigma_{zz}")

str = {'\sigma_{rr} < \sigma_{\theta\theta} < \sigma_{zz}'};
ht = text(0.3,0.03,str,'Interpreter','tex');
set(ht,'Rotation',57)
set(ht,'FontSize',15)

str = {'\sigma_{\theta\theta} < \sigma_{rr} < \sigma_{zz}'};
ht = text(0.08,0.85,str,'Interpreter','tex');
set(ht,'Rotation',0)
set(ht,'FontSize',15)

str = {'\sigma_{\theta\theta} < \sigma_{zz}  < \sigma_{rr}'};
ht = text(0.2,1.5,str,'Interpreter','tex');
set(ht,'Rotation',0)
set(ht,'FontSize',15)

str = {' \sigma_{zz} < \sigma_{\theta\theta} < \sigma_{rr}'};
ht = text(1.3,1.4,str,'Interpreter','tex');
set(ht,'Rotation',57)
set(ht,'FontSize',15)

str = {' \sigma_{zz} <  \sigma_{rr} < \sigma_{\theta\theta} '};
ht = text(1.2,1.15,str,'Interpreter','tex');
set(ht,'Rotation',0)
set(ht,'FontSize',15)

str = {'  \sigma_{rr} <  \sigma_{zz} < \sigma_{\theta\theta} '};
ht = text(1,0.5,str,'Interpreter','tex');
set(ht,'Rotation',0)
set(ht,'FontSize',15)

xlim([0,2])
ylim([0,2]);

%% Failure envelope
subplot(122)
plot(xl,[1,1],'Color', clrs(1,:)); hold on;
plot(xl,xl,'Color', clrs(2,:))
plot(xl,2*xl-1,'Color', clrs(3,:))
grid on;
plot_MC(1/2.7, 0, 0,[0 0 0])
xlim([0,3]);
ylim([0,3.5]);
xl = xlim;



function [h] = plot_MC(pp, taup, C_over_sigz,clr,symbl,lw)
% taup: tau/p, the shear tractions
% rR: r/R, non-deminsional conduit radius
% C_over_sigz: UCS divided by vertical stress
% symbl: symbol for plotting
% lw: line weight for plotting
% pp: pore pressure gradient (Pp/sigma_z)

phi = 35/180*pi;
c = C_over_sigz;
q = tan(pi/4 + 1/2*phi)^2;

% These expressions are derived analytically in the derivations.m script
% TODO: the expressions have additional components that might come into
% play depending on taup. Add those.
bottom = @(taup,pp,c,q)   -(c + pp - c*q - 2*pp*q + pp*q^2 + q*(4*c^2*taup^2 + c^2 - 8*c*pp*q*taup^2 - 2*c*pp*q + 8*c*pp*taup^2 + 2*c*pp + 4*c*q*taup^2 + 2*c*q - 4*c*taup^2 - 2*c + 4*pp^2*q^2*taup^2 + pp^2*q^2 - 8*pp^2*q*taup^2 - 2*pp^2*q + 4*pp^2*taup^2 + pp^2 - 4*pp*q^2*taup^2 - 2*pp*q^2 + 8*pp*q*taup^2 + 4*pp*q - 4*pp*taup^2 - 2*pp + q^2 - 4*q*taup^2 - 2*q + 1)^(1/2) - q^2 + (4*c^2*taup^2 + c^2 - 8*c*pp*q*taup^2 - 2*c*pp*q + 8*c*pp*taup^2 + 2*c*pp + 4*c*q*taup^2 + 2*c*q - 4*c*taup^2 - 2*c + 4*pp^2*q^2*taup^2 + pp^2*q^2 - 8*pp^2*q*taup^2 - 2*pp^2*q + 4*pp^2*taup^2 + pp^2 - 4*pp*q^2*taup^2 - 2*pp*q^2 + 8*pp*q*taup^2 + 4*pp*q - 4*pp*taup^2 - 2*pp + q^2 - 4*q*taup^2 - 2*q + 1)^(1/2) - 1)/(2*(q^2*taup^2 + 2*q*taup^2 + q + taup^2));
top = @(taup,pp,c,q) (c*q - pp - c + 2*pp*q - pp*q^2 + q*(4*c^2*taup^2 + c^2 - 8*c*pp*q*taup^2 - 2*c*pp*q + 8*c*pp*taup^2 + 2*c*pp + 4*c*q*taup^2 + 2*c*q - 4*c*taup^2 - 2*c + 4*pp^2*q^2*taup^2 + pp^2*q^2 - 8*pp^2*q*taup^2 - 2*pp^2*q + 4*pp^2*taup^2 + pp^2 - 4*pp*q^2*taup^2 - 2*pp*q^2 + 8*pp*q*taup^2 + 4*pp*q - 4*pp*taup^2 - 2*pp + q^2 - 4*q*taup^2 - 2*q + 1).^(1/2) + q^2 + (4*c^2*taup^2 + c^2 - 8*c*pp*q*taup^2 - 2*c*pp*q + 8*c*pp*taup^2 + 2*c*pp + 4*c*q*taup^2 + 2*c*q - 4*c*taup^2 - 2*c + 4*pp^2*q^2*taup^2 + pp^2*q^2 - 8*pp^2*q*taup^2 - 2*pp^2*q + 4*pp^2*taup^2 + pp^2 - 4*pp*q^2*taup^2 - 2*pp*q^2 + 8*pp*q*taup^2 + 4*pp*q - 4*pp*taup^2 - 2*pp + q^2 - 4*q*taup^2 - 2*q + 1).^(1/2) + 1)/(2*(q^2*taup^2 + 2*q*taup^2 + q + taup^2));
rightminus = @(k,taup,pp,c,q) -(2*c - 4*k + 2*pp + q + c*q - 2*k*q - pp*q - q*(4*c^2*taup^2 + c^2 - 16*c*k*taup^2 - 4*c*k - 8*c*pp*q*taup^2 - 2*c*pp*q + 8*c*pp*taup^2 + 2*c*pp + 4*c*q*taup^2 + 2*c*q + 2*c + 16*k.^2*taup^2 + 4*k.^2 + 16*k*pp*q*taup^2 + 4*k*pp*q - 16*k*pp*taup^2 - 4*k*pp - 8*k*q*taup^2 - 4*k*q - 4*k + 4*pp^2*q^2*taup^2 + pp^2*q^2 - 8*pp^2*q*taup^2 - 2*pp^2*q + 4*pp^2*taup^2 + pp^2 - 4*pp*q^2*taup^2 - 2*pp*q^2 + 4*pp*q*taup^2 + 2*pp + q^2 + 2*q + 1).^(1/2) - pp*q^2 + q^2)./(2*(- q^2*taup^2 + q + 1));
rightplus = @(k,taup,pp,c,q) -(2*c - 4*k + 2*pp + q + c*q - 2*k*q - pp*q + q*(4*c^2*taup^2 + c^2 - 16*c*k*taup^2 - 4*c*k - 8*c*pp*q*taup^2 - 2*c*pp*q + 8*c*pp*taup^2 + 2*c*pp + 4*c*q*taup^2 + 2*c*q + 2*c + 16*k.^2*taup^2 + 4*k.^2 + 16*k*pp*q*taup^2 + 4*k*pp*q - 16*k*pp*taup^2 - 4*k*pp - 8*k*q*taup^2 - 4*k*q - 4*k + 4*pp^2*q^2*taup^2 + pp^2*q^2 - 8*pp^2*q*taup^2 - 2*pp^2*q + 4*pp^2*taup^2 + pp^2 - 4*pp*q^2*taup^2 - 2*pp*q^2 + 4*pp*q*taup^2 + 2*pp + q^2 + 2*q + 1).^(1/2) - pp*q^2 + q^2)./(2*(- q^2*taup^2 + q + 1));
left = @(k,taup,pp,c,q) (c + pp - q + 2*c*q + 2*k*q + pp*q - (4*c^2*taup^2 + c^2 + 16*c*k*q*taup^2 + 4*c*k*q - 8*c*pp*q*taup^2 - 2*c*pp*q + 8*c*pp*taup^2 + 2*c*pp - 2*c*q - 4*c*taup^2 - 2*c + 16*k.^2*q^2*taup^2 + 4*k.^2*q^2 - 16*k*pp*q^2*taup^2 - 4*k*pp*q^2 + 16*k*pp*q*taup^2 + 4*k*pp*q - 4*k*q^2 - 8*k*q*taup^2 - 4*k*q + 4*pp^2*q^2*taup^2 + pp^2*q^2 - 8*pp^2*q*taup^2 - 2*pp^2*q + 4*pp^2*taup^2 + pp^2 + 2*pp*q^2 + 4*pp*q*taup^2 - 4*pp*taup^2 - 2*pp + q^2 + 2*q + 1).^(1/2) + 4*k*q^2 - 2*pp*q^2 - 1)./(2*(q^2 + q - taup^2));
tausingular = @(q) (q + 1)^(1/2)/q;
k = 0:.01:5;
kk = [k k]; % used for solutions with more than one k value

tsingular = tausingular(q);

topval = top(taup, pp, c, q);
bottomval = bottom(taup, pp, c, q);

right = [rightminus(k, taup, pp, c, q) rightplus(k, taup, pp, c, q) ];
left = left(k, taup, pp, c, q);

topline = ones(size(k))*topval;
bottomline = ones(size(k))*bottomval;

if taup > tsingular
    % need both solutions
    ki = k < eval(subs(kright));
    righti = [ki ki];
    righti = righti & (right > bottomval) & (right < topval);
    
else
    % only need the minus solution
    righti = [ones(size(k)) zeros(size(k))];
    righti = righti & (right > bottomval) & (right < topval);
    
end


lefti = (left > bottomline) & (left < topline);

toprightk = max(k(righti));
topleftk = max(k(lefti));
bottomrightk = min(k(righti));
bottomleftk = min(k(lefti));

topi = (k < toprightk) & (k > topleftk);
bottomi = (k < bottomrightk) & (k > bottomleftk);

%plot(k, subs(rightplus)); hold on
h = plot(kk(righti), right(righti), 'Color', clr, 'DisplayName', string(taup)); hold on;
plot(k(topi), topline(topi), 'Color', clr); hold on;
plot(k(bottomi), bottomline(bottomi),'Color',  clr);

%plot(k, subs(leftplus)); hold on
plot(k(lefti), left(lefti),'Color',  clr); hold on;
%xlim([0,4])
%ylim([0,4])
grid on
xlabel('k')
ylabel('p/\sigma_{zz}')

end

