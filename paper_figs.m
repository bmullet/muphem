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

clrs = parula(3);


set(0,'defaultFigurePosition', [defpos(1) defpos(2) 2.2*width*100, height*100]);

rR = 1;
lambda = linspace(0.01,2,500);
cond1 = @(lambda, rR) lambda;
cond2 = @(lambda, rR) lambda*(rR^2 + 1) - rR^2;
cond3 = @(lambda, rR) lambda*(1 - rR^2) + rR^2;

subplot(121)
plot(lambda, cond1(lambda, rR), 'Color', clrs(1,:)); hold on;
plot(lambda, cond2(lambda, rR), 'Color', clrs(2,:));
plot(lambda, cond3(lambda, rR), 'Color', clrs(3,:));

xlabel('S/\sigma_{zz}');
ylabel('p/\sigma_{zz}');
ylim([0,2])
xlim([0,2])
grid on

str = {'\sigma_{rr} < \sigma_{\theta\theta} < \sigma_{zz}'};
ht = text(0.3,0.03,str,'Interpreter','tex');
set(ht,'Rotation',57)
set(ht,'FontSize',15)

% str = {'\sigma_{\theta\theta} < \sigma_{rr} < \sigma_{zz}'};
% ht = text(0.1,0.5,str,'Interpreter','tex');
% set(ht,'Rotation',37)
% set(ht,'FontSize',15)

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

%% r/R = 1.1

subplot(122)
rR = 1;
plot(lambda, cond1(lambda, rR), 'Color', clrs(1,:)); hold on;
plot(lambda, cond2(lambda, rR), 'Color', clrs(2,:));
plot(lambda, cond3(lambda, rR), 'Color', clrs(3,:));

rR = 1.3;
lambda = linspace(0.01,2,500);
plot(lambda, cond1(lambda, rR), 'Color', clrs(1,:), 'LineStyle', '--','LineWidth', 3); hold on;
plot(lambda, cond2(lambda, rR), 'Color', clrs(2,:), 'LineStyle', '--','LineWidth', 3);
plot(lambda, cond3(lambda, rR), 'Color', clrs(3,:), 'LineStyle', '--','LineWidth', 3);

rR = 2;
lambda = linspace(0.01,2,500);
plot(lambda, cond1(lambda, rR), 'Color', clrs(1,:), 'LineStyle', ':'); hold on;
plot(lambda, cond2(lambda, rR), 'Color', clrs(2,:), 'LineStyle', ':');
plot(lambda, cond3(lambda, rR), 'Color', clrs(3,:), 'LineStyle', ':');

xlabel('S/\sigma_{zz}');
ylabel('p/\sigma_{zz}');
ylim([0,2])
xlim([0,2])
grid on

p1 = plot(xlim, [-1 -1], '-k');
p2 = plot(xlim, [-1 -1], '--k','LineWidth', 2);
p3 = plot(xlim, [-1 -1], ':k');
legend([p1, p2, p3],'r/R = 1', 'r/R = 1.3', 'r/R = 2')

%% Anotations
subplot(121); hold on;
str = {'(a)'};
ht = text(0.03,42,str,'Interpreter','tex');
set(ht,'Rotation',0)
set(ht,'FontSize',20)

subplot(122); hold on;
str = {'(b)'};
ht = text(0.03,42,str,'Interpreter','tex');
set(ht,'Rotation',0)
set(ht,'FontSize',20)

%% Failure
clrs = parula(3);
phi = 38/180*pi;
figure

set(0,'defaultFigurePosition', [defpos(1) defpos(2) 2.2*width*100, height*100]);

subplot(121)
rR = 1;
plot_MC(1,0,'-',3);
plot_MC(1.2,0,':',2);
p2 = plot(xlim, [-1 -1], ':k','LineWidth',2);
p1 = plot(xlim, [-1 -1], '-k','LineWidth',3);



C = 0;
%b = 2*sin(phi)/(sqrt(3)*(3-sin(phi))); a = 6*C*cos(phi)/(sqrt(3)*(3-sin(phi))); % circumscribe
b = sin(phi)/sqrt(9 + 3*sin(phi)^2);   a = 6*C*cos(phi)/(sqrt(3)*(3+sin(phi))); % inscribes
%b = 2*sin(phi)/(sqrt(3)*(3+sin(phi))); a = 3*C*cos(phi)/(sqrt(9+3*sin(phi))^2); % middle circumscribes


colors = {'k','b','g','r'};
lambda = 0:.0001:5;
d = [500000];
for j = 1:1
sigz = 2700*9.8*d(j);

fh = @(lam) (2*lam*sigz + (4*a^2 + 16*a*b*lam*sigz + 8*a*b*sigz + 16*b^2*lam^2*sigz^2 + 16*b^2*lam*sigz^2 + 4*b^2*sigz^2 - (4*lam^2*sigz^2)/3 + (8*lam*sigz^2)/3 - (4*sigz^2)/3)^(1/2))/(2*sigz);
fl = @(lam) (2*lam*sigz - (4*a^2 + 16*a*b*lam*sigz + 8*a*b*sigz + 16*b^2*lam^2*sigz^2 + 16*b^2*lam*sigz^2 + 4*b^2*sigz^2 - (4*lam^2*sigz^2)/3 + (8*lam*sigz^2)/3 - (4*sigz^2)/3)^(1/2))/(2*sigz);

pzp = nan(size(lambda));
pzm = nan(size(lambda));

for i=1:length(lambda)
    try
        pzp(i) =  fh(lambda(i));
        pzm(i) =  fl(lambda(i));
    catch
        1;
    end
end

pzp(imag(pzp) ~= 0) = nan;
pzm(imag(pzm) ~= 0) = nan;
pz = [fliplr(pzp), pzm];
lplot = [fliplr(lambda), lambda];

plot(lambda, pzp,strcat('-', colors{j})); hold on
 plot(lambda, pzm,strcat('-', colors{j}));
%plot(lplot, pz,strcat('-', colors{j})); hold on
% ylim([0,2])
% xlim([0,2])
end
%legend('\infty','5000','1000','500')

legend([p1,p2],'r/R = 1', 'r/R = 1.2','Location','Southeast')

subplot(122)
rR = 1;
plot_MC(1,0,'-',3);
plot_MC(1,0.5,':',2);
p2 = plot(xlim, [-1 -1], ':k','LineWidth',2);
p1 = plot(xlim, [-1 -1], '-k','LineWidth',3);
legend([p1,p2],'C/\sigma_{zz} = 0', 'C/\sigma_{zz} = 0.5','Location','Southeast')

subplot(121); hold on;
str = {'(a)'};
ht = text(0.15,4.7,str,'Interpreter','tex');
set(ht,'Rotation',0)
set(ht,'FontSize',20)

subplot(122); hold on;
str = {'(b)'};
ht = text(0.15,4.7,str,'Interpreter','tex');
set(ht,'Rotation',0)
set(ht,'FontSize',20)

%% Change of r-z failure with tau

set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*100]);


tau_flag = 'sigz'; % Can be 'p' or 'sigz' for how to normalize

clrs = parula(3);
s = 0.05:0.01:1;
fades = [0.25, 1]; % fade factor
taus = [0.0, 0.4];

if tau_flag == 'p'
 fail_rz = @(p,tau,c,q)   .5*(p+1+((p-1)^2 + 4*(tau*p)^2)^0.5) - c - q*0.5*(p+1-((p-1)^2 + 4*(tau*p)^2)^0.5);
 fail_rt = @(p,s,tau,c,q) - 2*s + p + c + q*0.5*(p + 1 - ((p-1)^2 + 4*(tau*p)^2)^.5);
 fail_tz = @(p,s,tau,c,q) 0.5*(p + 1 + ((p-1)^2 + 4*(tau*p)^2)^.5) - c - q*(2*s-p);
else
 fail_rz = @(p,tau,c,q)   .5*(p+1+((p-1)^2 + 4*(tau)^2)^0.5) - c - q*0.5*(p+1-((p-1)^2 + 4*(tau)^2)^0.5);
 fail_rt = @(p,s,tau,c,q) - 2*s + p + c + q*0.5*(p + 1 - ((p-1)^2 + 4*(tau)^2)^.5);
 fail_tz = @(p,s,tau,c,q) 0.5*(p + 1 + ((p-1)^2 + 4*(tau)^2)^.5) - c - q*(2*s-p);
end

c = 0;
phi = 38/180*pi;
%cohesion = 5e6; mu = tan(phi); C = 2*cohesion*((mu^2 + 1)^(1/2) + mu);

q = tan(pi/4 + 1/2*phi)^2;

for ii = 1:2

tau = taus(ii);
fade = fades(ii);
    
frz = ones(size(s));
frz_lower = frz * fzero(@(x) fail_rz(x,tau,c,q), 0.237);
frz_upper = frz * fzero(@(x) fail_rz(x,tau,c,q), 5);
frt = nan(size(s));
ftz = nan(size(s));

for i = 1:length(s)
    
    try
        frt(i) = fzero(@(x) fail_rt(x, s(i), tau, c, q), [-1, 1]);
        
    catch ME
        disp('no soln')
    end
    
    try
        ftz(i) = fzero(@(x) fail_tz(x, s(i), tau, c, q), [-1, 1]);
    catch ME
        disp('no soln')
    end
end

i1 = (frz_lower > frt) & ((frz_lower < ftz) + isnan(ftz));
i2 = (frz_upper > frt) & ((frz_upper < ftz) + isnan(ftz));

p1 = plot(s(i1), frz_lower(i1),':k'); hold on;
p1.Color(4) = fade;
p1 = plot(s(i2), frz_upper(i2),':k');
p1.Color(4) = fade;



i3 = (frt > frz_lower) & (frt < frz_upper);
i4 = (ftz > frz_lower) & (ftz < frz_upper);
p1 = plot(s(i3), frt(i3),'Color', [0 0 0], 'LineStyle', ':');
p1.Color(4) = fade;
p1 = plot(s(i4), ftz(i4),'Color', [0 0 0], 'LineStyle', ':');
p1.Color(4) = fade;
xlim([0,1])
ylim([0,1])
xlabel('S/\sigma_{zz}')
ylabel('p/\sigma_{zz}')

x = [s(i1), s(i3), fliplr(s(i2)), fliplr(s(i4))];
y = [frz_lower(i1), frt(i3), fliplr(frz_upper(i2)), fliplr(ftz(i4))];

if ii==2
    p1 = patch(x,y,'g','FaceAlpha',.1,'LineStyle', 'none');
%     p1.FaceVertexAlphaData = 0.01;    % Set constant transparency 
%     p1.FaceAlpha = 'flat' ;          % Interpolate to find face transparency
end

% with tau

%figure

if tau_flag == 'sigz'
%bound_rt = @(p,s,tau) p*(3*s-p-1) - (s*(2*s-1) - 1/2*(tau)^2);
bound_rt = @(p,s,tau) 2*s - p - (1/2*(1+p) - sqrt((1/2*(p-1))^2 + (tau)^2));
%bound_rz = @(p,s,tau) (4*s-3*p-1)^2 - ((p-1)^2 + 4*(tau)^2);
bound_rz = @(p,s,tau) 2*s - p - (1/2*(1+p) + sqrt((1/2*(p-1))^2 + (tau)^2));
else
   
bound_rt = @(p,s,tau) 2*s - p - (1/2*(1+p) - sqrt((1/2*(p-1))^2 + (tau*p)^2));
bound_rz = @(p,s,tau) 2*s - p - (1/2*(1+p) + sqrt((1/2*(p-1))^2 + (tau*p)^2));
end


prt = nan(size(s));
pzt = nan(size(s));

for i = 1:length(s)
   prt(i) = fzero(@(x) bound_rt(x,s(i),tau),s(i)+tau);
   pzt(i) = fzero(@(x) bound_rz(x,s(i),tau),2*s(i) - 1-tau);
end

p1 = plot(s,prt,'Color', clrs(1,:), 'LineStyle', '-'); hold on
p1.Color(4) = fade;
p1 = plot(s,pzt,'Color', clrs(2,:), 'LineStyle', '-');
p1.Color(4) = fade;
ylim([0,1])
xlim([0,1])

end

%%
set(0,'defaultFigurePosition', [defpos(1) defpos(2) 2.2*width*100, height*100]);

taus = [0.0, 0.01, 0.2, 1, 5, 50];

clrs = parula(length(taus));

for ii = 1:length(taus)
tp = taus(ii);
    pz = 0:.001:.99999;
th = nan(size(pz));

for i=1:length(pz)
    th(i) = 1/2*atand(2*tp*pz(i)/(pz(i)-1));
end

subplot(121)
plot(pz,abs(th),'Color',clrs(ii,:),'DisplayName',"\tau/p = " + num2str(tp)); hold on;
xlabel('p/\sigma_{zz}')
ylabel('Rotation (degrees)')


end
ylim([0,45])
legend('Location','Southeast')
    
%% t-star

tstar = @(c,q)  ((1./((c + q).*(1-c))).^(1/2).*(c + q - 1))/2;

phi = [0:.1:45]*pi/180;
dt = 0.001;
c = 0:.01:1-dt;
q = tan(pi/4 + 1/2*phi).^2;

[C,Q] = meshgrid(c,q);
[~,PHI] = meshgrid(c,phi);

ts = tstar(C,Q);
subplot(122)
[C,p1] = contour(C,PHI/pi*180,abs(ts), [ 0.1, 0.25, 0.5, 1, 2]);
p1.LineWidth = lw;
clabel(C,p1,'FontSize',17,'Color','black')
xlabel('C/\sigma_{zz}')
ylabel('\phi (degrees)')
% set(p1, 'EdgeColor','none')
% colorbar
% colormap('jet')
% caxis([0,2])

%semilogy(c,abs(ts)); hold on;
%xlabel('c/\sigma_{zz}')
%ylabel('\tau*')
%ylim([1e-1, 1e1])


%% Anotations
subplot(121); hold on;
str = {'(a)'};
ht = text(0.03,42,str,'Interpreter','tex');
set(ht,'Rotation',0)
set(ht,'FontSize',20)

subplot(122); hold on;
str = {'(b)'};
ht = text(0.03,42,str,'Interpreter','tex');
set(ht,'Rotation',0)
set(ht,'FontSize',20)

%% MSH plot
load('paperfigs/MSH.mat');

A = out{1}; zvec = out{2}; pvec = out{3}; ugvec = out{4}; umvec = out{5};
phivec = out{6}; rhogvec = out{7}; chidvec = out{8};

[Srr, Szz, Stt, Srz] = kirsch (zvec,pvec,A,ugvec,umvec,rhogvec,phivec,pvec);

A.mc.phi = 30/180*pi;

plotfailureprofiles(A,pvec,Szz,nan(size(Szz)),Srz,zvec,pvec,true);


%% Critical radius
noshear1 = importdata('CriticalRadius/vary_lambda_constant_p_noshear');
noshear2 = importdata('CriticalRadius/vary_lambda_constant_p_noshear_under7');
withshear = importdata('CriticalRadius/vary_lambda_constant_p_with_shear_5to1');
[~,idx1] = min(abs(noshear1(1,:)-0.658));
[~,idx2] = min(abs(noshear2(1,:)-0.658));


plot(noshear1(1,idx1:end),noshear1(2,idx1:end),'-r')
hold on
h2 = plot(noshear2(1,1:idx2),noshear2(2,1:idx2),'-r');

h1 = plot(withshear(1,:),withshear(2,:), '-b');
xlim([0.5,0.75]);

legend([h1,h2], 'With shear','No shear','Location','northwest');
xlabel('S/\sigma_{zz}');
ylabel('Min. stable radius (m)')

str = {'  \sigma_{rr} <  \sigma_{\theta\theta} < \sigma_{zz}  '};
ht = text(0.5,60,str,'Interpreter','tex');
set(ht,'Rotation',0)
set(ht,'FontSize',12)

str = {'    \sigma_{rr} < \sigma_{zz} < \sigma_{\theta\theta}  '};
ht = text(0.62,70,str,'Interpreter','tex');
set(ht,'Rotation',35)
set(ht,'FontSize',12)

str = {'  \sigma_{rr} <  \sigma_{\theta\theta} < \sigma_{zz}  '};
ht = text(0.52,120,str,'Interpreter','tex');
set(ht,'Rotation',0)
set(ht,'FontSize',12)

str = {'    \sigma_{rr} < \sigma_{zz} < \sigma_{\theta\theta}  '};
ht = text(0.65,125,str,'Interpreter','tex');
set(ht,'Rotation',40)
set(ht,'FontSize',12)


%% Failure depth
lw = 3;

plastic = importdata('MSHFigures/plastic_failure.txt');
elastic = importdata('MSHFigures/elastic_solution.txt');

phi = A.mc.phi;
cohesion = A.mc.C; mu = tan(phi);

C = 2*cohesion*((mu^2 + 1)^(1/2) + mu);
q = tan(pi/4 + 1/2*phi)^2;

rzrf = @(p, sigz, lambda, beta) sqrt((p./sigz - lambda)./(1/q - C./(q*sigz) - lambda));

rfailzr = rzrf(Srr, Szz, A.lambda);


hold on
rp = plastic(:,1);
zp = plastic(:,2);
[zp,idx] = sort(zp);
rp = rp(idx);

re = elastic(:,1);
ze = elastic(:,2);
[ze,idx] = sort(ze);
re = re(idx);


p2 = plot(re,ze,'-r');

p3 = plot(rp,zp,'-k');
p1 = plot(rfailzr*A.r, zvec,'--b','DisplayName', 'Analytical','LineWidth',3);

legend([p1,p2,p3],'Analytical','FEM Elastic','FEM Elastoplastic ')

xlim([40,50])
ylim([-3001, 0])

xlabel('r (m)'); ylabel('z (m)')



%%
% Eruption progression

clrs = parula(3);
phi = 38/180*pi;
figure

set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*100]);

rR = 1;
plot_MC(1,0,'-',3);

p2 = plot(xlim, [-1 -1], ':k','LineWidth',2);


C = 0;
b = sin(phi)/sqrt(9 + 3*sin(phi)^2);   a = 6*C*cos(phi)/(sqrt(3)*(3+sin(phi))); % inscribes

colors = {'k','b','g','r'};
lambda = 0:.0001:5;
d = [500000];

xlim([0,1])
ylim([0,1])


% define conduit regions
x = [0.4, 0.5, 0.5, 0.4];
y = [0.45, 0.45, 0.95, 0.95];

p1 = patch(x,y,'g','FaceAlpha',.4,'LineStyle', 'none');

x = [0.4, 0.5, 0.5, 0.4];
y = [0.05, 0.05, 0.45, 0.45];

p2 = patch(x,y,'b','FaceAlpha',.4,'LineStyle', 'none');

legend([p1,p2],'Early Widening','Late Collapse','Location','southeast','FontSize',15)

%%
% Stable radius
d1 = importdata('CriticalRadius2/vary_lambda_constant_p_no_shear_5to1');
d2 = importdata('CriticalRadius2/vary_lambda_constant_phi65_2');
d3 = importdata('CriticalRadius2/vary_lambda_constant_phi67_2');
d4 = importdata('CriticalRadius2/vary_lambda_constant_ezz_5to1');

plot(d3(1,:),d3(2,:)); hold on;
plot(d2(1,:),d2(2,:));
plot(d4(1,:),d4(2,:));
plot(d1(1,:),d1(2,:));

legend('\phi_f = 0.67', '\phi_f = 0.65', 'shear strain','no shear')

xlabel('S/\sigma_{zz}')
ylabel('Min. stable radius (m)')

function plot_MC(rR,C_over_sigz,symbl,lw)
clrs = parula(3);
phi = 38/180*pi;
cohesion = 5e6; mu = tan(phi); C = 2*cohesion*((mu^2 + 1)^(1/2) + mu);

q = tan(pi/4 + 1/2*phi)^2;
sigz = 2700*9.8*500;
C = C_over_sigz*sigz;


lambda = linspace(0.01,5,500);
bound_1 = @(lambda, rR) lambda;
bound_2 = @(lambda, rR) lambda*(rR^2 + 1) - rR^2;
bound_3 = @(lambda, rR) lambda*(1 - rR^2) + rR^2;

A_failure = @(x,c,sig_z,rR) rR^2*(1/q-c/(q*sig_z)) + x*(1-rR^2);
B_failure = @(x,c,sig_z,rR) x/(q+1)*(rR^2+1) - c/sig_z * 1/(q+1)*rR^2 - q/(q+1)*x*(rR^2-1);
C_failure = @(x,c,sig_z,rR) x*(rR^2 + 1) - rR^2*(c/sig_z + q);
D_failure = @(x,c,sig_z,rR) x*(1-rR^2) + rR^2*(c/sig_z + q);
E_failure = @(x,c,sig_z,rR) 1/(1+q)*(c/(sig_z)*rR^2 + x*(q*(rR^2 + 1) + (1-rR^2)));
F_failure = @(x,c,sig_z,rR) rR^2*(c/(q*sig_z)- 1/q) + x*(1+rR^2) ; 

x1 = fzero(@(x) A_failure(x,C,sigz,rR) - bound_1(x, rR), 0.6);
x2 = fzero(@(x) B_failure(x,C,sigz,rR) - bound_2(x, rR), 0.6);
x3 = fzero(@(x) C_failure(x,C,sigz,rR) - bound_3(x, rR), 0.6);
x4 = fzero(@(x) D_failure(x,C,sigz,rR) - bound_1(x, rR), 0.6);
x5 = fzero(@(x) E_failure(x,C,sigz,rR) - bound_2(x, rR), 0.6);
x6 = fzero(@(x) F_failure(x,C,sigz,rR) - bound_3(x, rR), 0.6);

plot(lambda, bound_1(lambda, rR), 'Color', clrs(1,:),'LineWidth',lw, 'LineStyle', symbl); hold on;
plot(lambda, bound_2(lambda, rR), 'Color', clrs(2,:),'LineWidth',lw, 'LineStyle', symbl);
plot(lambda, bound_3(lambda, rR), 'Color', clrs(3,:),'LineWidth',lw, 'LineStyle', symbl);


%plot([no_cohesion, no_cohesion], ylim, '--r')
plot([x1, x2], [A_failure(x1,C,sigz,rR), A_failure(x2,C,sigz,rR)], strcat(symbl,'k'),'LineWidth', lw)
plot([x2, x3], [B_failure(x2,C,sigz,rR), B_failure(x3,C,sigz,rR)], strcat(symbl,'k'),'LineWidth', lw)
plot([x3, x4], [C_failure(x3,C,sigz,rR), C_failure(x4,C,sigz,rR)], strcat(symbl,'k'),'LineWidth', lw)
plot([x4, x5], [D_failure(x4,C,sigz,rR), D_failure(x5,C,sigz,rR)], strcat(symbl,'k'),'LineWidth', lw)
plot([x5, x6], [E_failure(x5,C,sigz,rR), E_failure(x6,C,sigz,rR)], strcat(symbl,'k'),'LineWidth', lw)
plot([x1, x6], [F_failure(x1,C,sigz,rR), F_failure(x6,C,sigz,rR)], strcat(symbl,'k'),'LineWidth', lw)
xlabel('S/\sigma_{zz}');
ylabel('p/\sigma_{zz}');
ylim([0,5])
xlim([0,5])
end
