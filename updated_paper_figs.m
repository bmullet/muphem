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

%%
% Get data for comsol
%load("example_eruption.mat")
% srz and p
mat = [out{2} out{17}];
writematrix(mat, 'p50.csv')
mat = [out{2} out{18}];
writematrix(mat, 'srz50.csv')

% cohesion
mat = [out{2} A.mc.C(out{2})];
writematrix(mat, 'cohesion50.csv')

%% Comparison to COMSOL
lw = 3;

load("comsol_failure.mat")

Srr = out{17};
Srz = out{18};
zvec = out{2};
A = out{1};
Szz = A.k.rho*9.8*abs(zvec);

with_tau = importdata('../COMSOL/mc_failure_with_tau.txt');
z_tau = with_tau(:,2);
r_tau = with_tau(:,1);
no_tau = importdata('../COMSOL/mc_failure_no_tau.txt');
z_no_tau = no_tau(:,2);
r_no_tau = no_tau(:,1);

epe = importdata('../COMSOL/epe0.txt');
z_epe = epe(:,2);
r_epe = epe(:,1);


phi = 30/180*pi;
pp = abs(zvec)*9.8*1000;
cohesion = A.mc.C(zvec); 
mu = tan(phi);

C = 2*cohesion.*((mu^2 + 1)^(1/2) + mu);

q = tan(pi/4 + 1/2*phi)^2;

%rzrf = @(p, sigz, C, lambda) sqrt((p./sigz - lambda)./(1/q - C./(q*sigz) - lambda));
rzrf = @(p, sigz, C, k, pp) (-(q.*(p - k.*sigz))./(C + pp - sigz - pp.*q + k.*q.*sigz)).^(1/2);

rfailzr = rzrf(Srr, Szz, C, A.lambda(zvec),pp);

hold on
[z_tau,idx] = sort(z_tau);
r_tau = r_tau(idx);

[z_no_tau,idx] = sort(z_no_tau);
r_no_tau = r_no_tau(idx);

%subplot(121)
% radial comparison of failure
r_epe = r_epe(z_epe>-1800);
z_epe = z_epe(z_epe>-1800);

[z_epe,idx] = sort(z_epe);
r_epe = r_epe(idx);
rr = smooth(r_epe, 50);
%p3 = plot(r_epe,z_epe,'-k','LineWidth',1); hold on;
p3 = plot(rr,z_epe,'-k','LineWidth',1); hold on;
p1 = plot(rfailzr*A.r, zvec,'--r','DisplayName', 'Analytical','LineWidth',2);

legend([p1,p3],'Analytical','FEM')

xlim([50,60])
ylim([-2001, -1000])

xlabel('r (m)'); ylabel('z (m)')

grid on;

% subplot(122)
% % von mises stress
% vm = importdata('../COMSOL/von_mises.txt');
% sigma_theta = 2*Szz.*A.lambda(zvec) - Srr;
% p = 1/3*(sigma_theta + Szz + Srr);
% von_mises = (3/2*((sigma_theta -p).^2 + (Srr - p).^2 + (Szz - p).^2 + 2*(Srz).^2)).^(1/2);
% plot(zvec, von_mises); hold on;
% plot(vm(:,1), vm(:,2))

%%
dstress = importdata("../COMSOL/diff_stress.txt");
plot(zvec, Srr - Szz); hold on;
plot(dstress(:,1), dstress(:,2));

%%
MC = Szz - pp - (C + q.*(Srr - pp)); hold on;
plot(MC, zvec);
plot([0, 0], ylim, '--r')

%% Example eruption
load("example_eruption.mat")
colors = lines(4);
colors = colors(3:end, :);
A = out{1};
zvec = out{2};
pvec = out{3};
ugvec = out{4};
umvec = out{5};
phivec = out{6};
Srr = out{17};
Srz = out{18};

figure
subplot(1,3,1);
plot(pvec, zvec, "Color", colors(1,:)); hold on;
plot(Srz, zvec, "Color", colors(2,:)); hold on;
plot(xlim, [-2252,-2252], '--r', "Linewidth", 1)
legend("p", "\tau");
grid on;
xlabel("Pa")
ylabel("z (m)")
subplot(1,3,3);
semilogx(ugvec, zvec, "Color", colors(1,:)); hold on;
semilogx(umvec, zvec, "Color", colors(2,:)); hold on;
grid on;
plot(xlim, [-2252,-2252], '--r', "Linewidth", 1)
legend("u_g", "u_m");
xlabel("m/s");
ylabel("z (m)");
set(gca, 'XTick', [10, 100])
subplot(132);
plot(Srz./Srr, zvec, "Color", colors(1,:)); hold on; grid on;
plot(xlim, [-2252,-2252], '--r', "Linewidth", 1)
xlabel("\tau/p");
ylabel("z (m)");



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

%% Minimum stable radius - vary k
figure
colors = parula(2)*0.9;

% r_min v. p
subplot(2,3,1);
load("phi_75_with_shear_constant_k.mat")
plot(lambdas, rvec, 'Color', colors(1,:)); grid on;
xlabel("p_{chamber}/\sigma_{zz}")
ylabel("R_{stable}")
% r_min v. k
subplot(234)
load("phi_75_with_shear_minimum_constant_p.mat")
plot(lambdas, rvec, 'Color', colors(1,:));
grid on;
xlabel("k")
ylabel("R_{stable}")

colors = parula(2)*0.9;
load("phi_75_with_shear_minimum.mat");
i = lambdas < 1.28;
%i(4:2:find(i,1,'last')-1) = 0;
subplot(2,3,[2,3,5,6])
plot(lambdas(i), rvec(i), 'Color', colors(1,:)); hold on;

load("phi_75_without_shear_minimum.mat");
plot(lambdas, rvec, '--', 'Color', colors(1,:), "Linewidth", 2);

load("phi_65_with_shear_minimum.mat");
plot(lambdas, rvec, 'Color', colors(2,:));

xlabel('k')
ylabel("R_{stable}")
xlim([0.5,1.4])

%grid on


[hleg,icons,plots] = legend("0.75", "0.75", "0.65");
title(hleg,'\phi_{frag}')
hleg.Title.Visible = 'on';
% the addition in height needed for the title:
title_hight = hleg.Position(4)/numel(plots);
hleg.Position([2 4]) = [hleg.Position(2)-title_hight hleg.Position(4)+title_hight];
% calculate new position for the elements in the legeng:
new_pos = fliplr(0.5/(numel(plots)+1):1/(numel(plots)+1):1);
hleg.Title.NodeChildren.Position = [0.5 new_pos(1) 0];
% set the text to the right position:
leg_txt = findobj(icons,'Type','Text');
txt_pos = cell2mat({leg_txt.Position}.');
txt_pos(:,2) = new_pos(2:end);
set(leg_txt,{'Position'},mat2cell(txt_pos,ones(numel(plots),1),3));
% set the icons to the right position:
leg_att = findobj(icons,'Type','Line');
% set(leg_att,{'YData'},mat2cell(repmat(repelem(new_pos(2:end).',...
%     numel(plots)),1,2),ones(numel(plots)*2,1),2))
set(leg_att,{'YData'},mat2cell(repmat(repelem(new_pos(2:end).',...
     2),1,2),ones(numel(plots)*2,1),2))

 
%% Stable radius eruption progression
load("phi_75_with_shear_constant_k.mat")
figure
subplot(1,2,1)
plot(lambdas, rvec)
grid on;
xlim([0.5, 0.8])
ylabel('R_{min}'); xlabel("p_{ch}/\sigma_{zz}")
legend("R_{min}")
subplot(1,2,2)
t = [0, 0.5, 1];
r = [0, 1, 1];
pch = [1,  0.7, 0.62];
plot(t, r, '--x','Color','k'); hold on;
plot(t, pch, '--x','Color','r');
xlabel("t/t_{collapse}")
ylabel("Parameter Value")
legend("R/R_{final}", "P_{chamber}/\sigma_{zz}")
ylim([0,1.1])


%% 

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

