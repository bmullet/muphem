%% Test 1 - Test for uniqueness
A = initA();
z = 1:A.depth;
plith = z*A.k.rho*A.g;
A.Pchamber = max(plith)-60e6;
p = [];
%vvec = [.01:.01:3 3:.1:11 11:.002:12];
vvec = [.01:.01:.2 .2:.1:10];
%vvec = [.001:.002:.05 .2:.1:10];

z = [];
for i = vvec
    A.v_chamber_i = i;
    [zvec,pvec,ugvec,umvec,phivec,rhogvec,chidvec,Qmvec,Qgvec,A] = incoodes(A);
    p = [p min(pvec)];
    z = [z max(zvec)];
    disp(i)
end
plot(vvec,p)
hold on
xlabel('v')
ylabel('p top')
plot(xlim,[1.01e5 1.01e5],'--r')

figure
plot(vvec,z)
xlabel('v')
ylabel('top z')

%% Test 2 - Vary k
k = logspace(-10,-13,10);
vi = [];
pt = [];
for i = k
    disp(['k = '])
    disp(i)
    out = muphem('multiflow2',10e6,i);
    p = out{3};
    vg = out{4};
    vi = [vi min(vg)];
    pt = [pt min(p)];
end
figure
plot(k,pt)
xlabel('k')
ylabel('pressure')
set(gca,'XScale','log')
figure
plot(k,vi)
xlabel('k')
ylabel('chamber v')
set(gca,'XScale','log')

%% Test permeability law
kw0 = 1e-13;
z = 0:10000;
A = initA();
ppore = A.rhow*A.g*z+A.Patm_; % hyrdostatic pore pressure
kw1 = A.kw0*exp(-ppore/A.Pstar);

kw2 = 10.^(-14 - 3.2*log10(z/1000+.2));
plot(kw1,z,kw2,z)
legend('k1','k2')
set(gca,'XScale','log','Ydir','reverse')

{A,zvec,pvec,ugvec,umvec,phivec,rhogvec,chidvec,Qmvec,Qgvec,failure,Sprincipal,slip};

%% Distribute solution
A = out{i}{1}; zvec = out{i}{2}; pvec = out{i}{3}; ugvec = out{i}{4}; umvec = out{i}{5};
phivec = out{i}{6}; rhogvec = out{i}{7}; chidvec = out{i}{8}; Qmvec = out{i}{9}; Qgvec = out{i}{10}; failure = out{i}{11}; Sprincipal = out{i}{12};
slip = out{i}{13};
disp('Output distributed')

%% Test 3 - Test for uniqueness
A = initA();
z = 1:A.depth;
plith = z*A.k.rho*A.g;
vvec = [.05:.2:15];
op = [-120e6 -100e6 -80e6 -40e6  -20e6 0];

p = zeros(length(op), length(vvec));
z = zeros(length(op), length(vvec));
for j = 1:length(op)
    A.Pchamber = max(plith)+op(j);
    for i = 1:length(vvec)
        A.v_chamber_i = vvec(i);
        try
            [zvec,pvec,ugvec,umvec,phivec,rhogvec,chidvec,Qmvec,Qgvec,A] = incoodes(A);
            p(j,i) = min(pvec);
            z(j,i) = max(zvec);
        catch
            p(j,i) = 10e6;
            z(j,i) = 500;
        end
        disp(vvec(i))
    end
    figure
    plot(vvec,p(j,:))
    hold on
    xlabel('v')
    ylabel('p top')
    plot(xlim,[1.01e5 1.01e5],'--r')
    title(['Overpressure = ' num2str(op(j)/1e6) ' MPa'])
    figure
    plot(vvec,z(j,:))
    xlabel('v')
    ylabel('top z')
    title(['Overpressure = ' num2str(op(j)/1e6) ' MPa'])
    drawnow()

end

%% Plot results of test three

figure
for i = 1:length(op)
    plot(vvec,p(i,:))
    hold on
end
xlabel('v')
ylabel('p top')
plot(xlim,[1.01e5 1.01e5],'--r')
legend('-120Mpa', '-100', '-80', '-40', '-20e6', '0')
title('Pressure at surface vs. initial velocity')


%% Test for interphase transition
ed = A.exdepth;
fd = A.fragdepth;
A.delF = 1;
bgi = (zvec>=ed&zvec<fd);
A.delF = 1;
delF = 1;
RHS1 = zeros(nnz(bgi),3);
phibg = phivec(bgi);
umbg = umvec(bgi);
ugbg = ugvec(bgi);
rhogbg = rhogvec(bgi);
pbg = pvec(bgi);
Qmbg = Qmvec(bgi); chidbg = chidvec(bgi);
zbg = zvec(bgi);
Fmgex = zeros(nnz(bgi),1);
F1 = Fmgex;
rbvec = Fmgex;
phitestex = Fmgex;
Is = Fmgex;
Fw = Fmgex;
Fmg1ex = Fmgex; Fmg2ex = Fmgex;
for i = 1:nnz(bgi)
    [Fmgex(i), phitestex(i), Fmg1ex(i), Fmg2ex(i)] = testfmg(zbg(i),[pbg(i), ugbg(i), umbg(i)],A);
end
fi = (zvec>=fd);
A.delF = 0;
delF = 0;
phif = phivec(fi);
umf = umvec(fi);
ugf = ugvec(fi);
rhogf = rhogvec(fi);
pf = pvec(fi);
Qmf = Qmvec(fi); chidf = chidvec(fi);
Fmgfr = zeros(nnz(fi),1);
zf = zvec(fi);
Fmgfr = zeros(nnz(fi),1);
F1 = Fmgfr;
Is = Fmgfr;
Fw = Fmgfr;
phitestfr = Fmgfr;
Fmg1fr = Fmgfr;
Fmg2fr = Fmgfr;
j = nnz(bgi);
for i=1:nnz(fi)
    [Fmgfr(i), phitestfr(i), Fmg1fr(i), Fmg2fr(i)] = testfmg(zf(i),[pf(i), ugf(i), umf(i)],A);
end
Fmg = [Fmgex; Fmgfr];
phitest = [phitestex; phitestfr];
Fmg1 = [Fmg1ex; Fmg1fr];
Fmg2 = [Fmg2ex; Fmg2fr];
testi = bgi | fi;

%% Bring down overpressure 

muvvec = [1e3, 1e4, 1e6, 1e8];
rvvec = [5, 25, 50, 100];
pvvec = [60e6:-30e6:-60e6];

Nmu = length(muvvec);
Nr = length(rvvec);
Np = length(pvvec);

it = 0;
N = Nmu*Nr*Np;   

vbot = nan(Nmu,Nr,Np);
vgtop = nan(Nmu,Nr,Np);
vmtop = nan(Nmu,Nr,Np);
ptop = nan(Nmu,Nr,Np);

for i = 1:length(muvvec)
    for j = 1:length(rvvec)
        for k = 1:length(pvvec)
            it = it+1;
            
            A.mu = muvvec(i);
            A.r = rvvec(j);
            p = pvvec(k);
            
            disp([num2str(it) ' out of ' num2str(N)]);
            disp(['Progress = ' num2str(it/N*100)]);
            disp(['mu = ' num2str(A.mu)]);
            disp(['p_ch = ' num2str(p)]);
            disp(['R = ' num2str(A.r)]);
            try
                out = muphem('multiflow2',p,A);
                A = out{1}; zvec = out{2}; pvec = out{3}; ugvec = out{4}; umvec = out{5};
                phivec = out{6}; rhogvec = out{7}; chidvec = out{8}; Qmvec = out{9}; Qgvec = out{10}; failure = out{11};
                disp('Output distributed')
                vbot(i,j,k) = min(umvec);
                vgtop(i,j,k) = max(ugvec);
                vmtop(i,j,k) = max(umvec);
                ptop(i,j,k) = min(pvec);
            catch
                disp('Code failed!')
            end
        end
    end
end
%% Make plots of various viscosities and failure points

A = initA();
muvvec = [1e3, 1e4, 1e6, 1e8];
r = 75;
pvvec = [60e6:-2e6:-60e6];
A.r = r;

Nmu = length(muvvec);
Np = length(pvvec);

it = 0;
N = Nmu*Np;   

vbot = nan(Nmu,Np);
vgtop = nan(Nmu,Np);
vmtop = nan(Nmu,Np);
ptop = nan(Nmu,Np);

for i = 1:length(muvvec)
        for k = 1:length(pvvec)
            it = it+1;
            
            A.mu = muvvec(i);
            p = pvvec(k);
            
            disp([num2str(it) ' out of ' num2str(N)]);
            disp(['Progress = ' num2str(it/N*100)]);
            disp(['mu = ' num2str(A.mu)]);
            disp(['p_ch = ' num2str(p)]);
            disp(['R = ' num2str(A.r)]);
            try
                out = muphem('multiflow2',p,A);
                A = out{1}; zvec = out{2}; pvec = out{3}; ugvec = out{4}; umvec = out{5};
                phivec = out{6}; rhogvec = out{7}; chidvec = out{8}; Qmvec = out{9}; Qgvec = out{10}; failure = out{11};
                disp('Output distributed')
                vbot(i,j,k) = min(umvec);
                vgtop(i,j,k) = max(ugvec);
                vmtop(i,j,k) = max(umvec);
                ptop(i,j,k) = min(pvec);
            catch
                disp('Code failed!')
            end
        end
end
save('failureandviscosity.mat')

%% Test for poissons ratio
A = initA();
for p = .1:.1:.5
    A.k.p = p;
    A.k.K = A.k.p/(1-A.k.p);
    out = muphem('multiflow2',0e6,A);
end


%% Plot results
%plotmuphem(A,zvec,pvec,ugvec,umvec,phivec,rhogvec,chidvec)

% Test for failure
[Srr, Szz, Stt] = kirsch(zvec,pvec,A);
[Smax,Sfail,failure] = mcfailure(A,Srr,Szz,Stt,zvec);
if (failure)
    disp('we have a failure!')
else
    disp('no failure')
end
plot_failure(zvec,pvec,phivec,A,Srr,Szz,Stt,Smax,Sfail,op)
disp('Pressure at conduit exit:')
disp(min(pvec));
vargout = {A,zvec,pvec,ugvec,umvec,phivec,rhogvec,chidvec,Qmvec,Qgvec,failure};

%% Plot stability tests
load('findstablerhgfunc.mat')
plot(hgvec,stabler(2,:),'Color',[148/255, 20/255, 155/255]);
hold on;
plot(hgvec,stabler(3,:),'Color',[32/255, 89/255, 204/255]);
plot(hgvec,stabler(4,:),'Color',[32/255, 204/255, 49/255]);
title('Stable Radius v. Water Content')
xlabel('Initial Water Mass Fraction')
ylabel('Minimum Stable Radius (m)')
legend('10^5','10^6','10^7')
set(gca,'yscale','log')
xlim([0.03,0.07])

figure
load('GOLD2.mat')
plot(muvec,stabler(:,1));
hold on;
plot(muvec,stabler(:,2));
plot(muvec,stabler(:,3));
title('Stable Radius v. Water Content')
xlabel('Initial Water Mass Fraction')
ylabel('Minimum Stable Radius (m)')
legend('0.03','0.05','0.07')
set(gca,'yscale','log')
set(gca,'xscale','log')

%% Test for failure
[Srr, Szz, Stt] = kirsch(zvec,pvec,A);
[Smax,Sfail,failure] = mcfailure(A,Srr,Szz,Stt,zvec);
if (failure)
    disp('we have a failure!')
else
    disp('no failure')
end
plot_failure(zvec,pvec,phivec,A,Srr,Szz,Stt,Smax,Sfail,op)
disp('Pressure at conduit exit:')
disp(min(pvec));
vargout = {A,zvec,pvec,ugvec,umvec,phivec,rhogvec,chidvec,Qmvec,Qgvec,failure};

%% Show plots of failure with changing overpressure
close all;
for i = 1:length(outvec)
    out = outvec{i};
    A = out{1}; zvec = out{2}; pvec = out{3}; ugvec = out{4}; umvec = out{5};
phivec = out{6}; rhogvec = out{7}; chidvec = out{8}; Qmvec = out{9}; Qgvec = out{10}; failure = out{11};
A.mu = @(x,y) A.mu;
disp('Output distributed')
% Test for failure
[Srr, Szz, Stt, Srz] = kirsch(zvec,pvec,A,ugvec,umvec,rhogvec,phivec,pvec);
[Smax,Sfail,failure] = mcfailure(A,Srr,Szz,Stt,Srz,zvec);
if (failure)
    disp('we have a failure!')
else
    disp('no failure')
end
plot_failure(zvec,pvec,phivec,A,Srr,Szz,Stt,Srz, Smax,Sfail)
pause

end

%% Plot exit velocities and viscosities and pressure
load('failureandviscosity3.mat')
subplot(1,2,1)
h = pcolor(pvvec./1e6,log10(muvvec(5:end)),vmtop(5:end,:));
set(h, 'EdgeColor', 'none');
ylabel('Log_{10} \mu (Pa-s)')
xlabel('Chamber overpressure (MPa)')
title('Melt exit velocity (m/s)')
colorbar
subplot(1,2,2)
h = pcolor(pvvec./1e6,log10(muvvec(5:end)),vgtop(5:end,:));
set(h, 'EdgeColor', 'none');
colorbar
xlabel('Chamber overpressure (MPa)')
ylabel('Log_{10} \mu (Pa-s)')
title('Gas exit velocity (m/s)')

%% Test exit velocities
v = .1:.1:15;
resid = zeros(size(v));
for i = 1:length(v)
    vi = v(i)
    resid(i) = matchPatm(vi, A);
end

%% Distribute solution
A = out{1}; zvec = out{2}; pvec = out{3}; ugvec = out{4}; umvec = out{5};
phivec = out{6}; rhogvec = out{7}; chidvec = out{8}; Qmvec = out{9}; Qgvec = out{10}; failure = out{11};
Srz = nan(size(phivec));
    mu = A.mu(phivec,pvec);
    Srz(zvec<A.fragdepth) = 4*mu(zvec<A.fragdepth).*umvec(zvec<A.fragdepth)/A.r;
    Srz(zvec>=A.fragdepth) = A.f0*rhogvec(zvec>=A.fragdepth).*ugvec(zvec>=A.fragdepth).^2./2;

%% Do MSH_Degruyter test
semilogx(pvz(:,1),pvz(:,2),'or')
set(gca,'YDir','Reverse')
ylabel('Depth (z)')
xlabel('Pressure (Pa)')
hold on;
semilogx(pvec,-zvec,'-b')
legend('From Degruyter (2012)', 'My implementation')

%% Quick test of stability
set(0,'defaultFigurePosition', [defpos(1)*2 defpos(2) width*100, height*100]);

lam=0.50; % lambda

sigz = -A.rhom0*zvec*9.8;
pz = pvec./sigz;

semilogx(pz, zvec);
xlim([1e-1, 1e1])
hold on
plot([0.23, 0.23],ylim, '--r')


%% Plot some failure things


phi = 38/180*pi;
cohesion = 5e6; mu = tan(phi); C = 2*cohesion*((mu^2 + 1)^(1/2) + mu);

q = tan(pi/4 + 1/2*phi)^2;
sigz = 2700*9.8*50000;

q1000 = 1/q*(1-C/(sigz));

lambda = linspace(0.01,10,500);
up_bound = @(lambda) lambda;
low_bound = @(lambda) 2*lambda-1;

right_failure = @(x,c,sig_z) 2/(q+1)*x - c/sig_z * 1/(q+1);
left_failure = @(x,c, sig_z) c/(q*sig_z) + 2*x - 1/q; 

xr = fzero(@(x) right_failure(x,C,sigz) - low_bound(x), 0.6);
xl = fzero(@(x) left_failure(x,C,sigz) - up_bound(x), 0.6);

no_cohesion = q1000;

plot(lambda, up_bound(lambda)); hold on;
plot(lambda, low_bound(lambda));
plot(xlim, [no_cohesion, no_cohesion], '--r')
plot([no_cohesion, no_cohesion], ylim, '--r')

plot([xr, 1], [right_failure(xr,C,sigz), right_failure(1,C,sigz)], '--g')
plot([xl, 1], [left_failure(xl,C,sigz), left_failure(1,C,sigz)], '--b')
xlabel('\lambda');
ylabel('p/\sigma_z');
% ylim([0,2])
% xlim([0,2])


%% dp failure TRY 2

%b = 2*sin(phi)/(sqrt(3)*(3-sin(phi))); a = 6*C*cos(phi)/(sqrt(3)*(3-sin(phi))); % circumscribe
b = sin(phi)/sqrt(9 + 3*sin(phi)^2);   a = 6*C*cos(phi)/(sqrt(3)*(3+sin(phi))); % inscribes
%b = 2*sin(phi)/(sqrt(3)*(3+sin(phi))); a = 3*C*cos(phi)/(sqrt(9+3*sin(phi))^2); % middle circumscribes


colors = {'k','b','g','r'};
d = [500000];
for j = 1:4
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
% plot(lambda, pzp,strcat('-', colors{j})); hold on
% plot(lambda, pzm,strcat('-', colors{j}));
plot(lplot, pz,strcat('-', colors{j})); hold on
% ylim([0,2])
% xlim([0,2])
end
legend('\infty','5000','1000','500')
xlabel('\lambda');
ylabel('p/\sigma_z');

%% drucker prager analytical
syms sigz p lmda
A = diag([sigz, p, 2*lmda*sigz - p]);
t = trace(A)/3;
d = A - t;
simplify(sqrt(sum(d.*d,'all'))/sigz)
simplify(sqrt(1/6*(expand((sigz - p)^2 + (p-(2*lmda*sigz-p))^2 + (sigz-(2*lmda*sigz-p))^2))))


%% General Failure Conditions
lcrit = 1/q*(1-C/sigz);
colors = {'k','b','r'};
draw = {'-','--','.'};
l = {.8, .3, lcrit};

rc = 50;
zrf = @(r, lambda) r.^2/rc^2*(1/q - C/(q*sigz)) + lambda*(1-r.^2/rc^2);
ztf = @(r, lambda) r.^2/rc^2*(C/(q*sigz) - 1/q) + lambda*(1+r.^2/rc^2);
trf = @(r, lambda) -r.^2/rc^2*(C/((q+1)*sigz)) + lambda/(1+q)*(1+r.^2/rc^2)...
    - q/(q+1)*lambda*(r.^2/rc^2-1);

R = linspace(50,100,100);
ax = [];
for i = 1:length(colors)
    ax(i) = plot(R,zrf(R,l{i}),strcat(draw{i}, colors{1})); hold on;
    plot(R,ztf(R,l{i}),strcat(draw{i}, colors{2})); hold on;
    plot(R,trf(R,l{i}),strcat(draw{i}, colors{3})); hold on;
end

ylim([0,1])
xlabel('r'); ylabel('p/\sigma_z');

legend(ax,string([l{:}]));
% dim = [0.2 0.5 0.3 0.3];
% str = {'Straight Line Plot','from 1 to 10'};
% annotation('textbox',dim,'String',str,'FitBoxToText','on','color','red');
words={'z/r', 'z/\theta', '\theta/r'};
for i = 1:length(color)  
    text(55, 1.01-(i*.07),words{i},'Color',colors{i},'FontSize',20);    
end

%% Failure distance

% A.rhom0 = 2600;
A.r = 100;

data = importdata('plam5.txt');
zvec = data(:,1);
pvec = data(:,2);

phi = 25/180*pi;
cohesion = 5e5; mu = tan(phi); C = 2*cohesion*((mu^2 + 1)^(1/2) + mu);

q = tan(pi/4 + 1/2*phi)^2;

%COMSOL
% C = 2.0527e6;
% q = 4.21346;

rzrf = @(p, sigz, lambda, beta) sqrt((1/beta*p./sigz - lambda)./(1/q - C./(q*sigz) - lambda));
rztf = @(p, sigz, lambda, beta) (1/beta*p./sigz - lambda)./(1/q - C./(q*sigz) + lambda);
lith = -A.rhom0*zvec*9.8;


beta = 1;
    rfailzr = rzrf(pvec, lith, 0.5, beta);
    plot(rfailzr*A.r, zvec, 'DisplayName', 'Analytical','LineWidth',3);
    hold on


% xlim([100, 115]);
% ylim([-2900, -1100]);
% xlim([80, 180])
% ylim([-3000, 0])
legend show


data = importdata('../Documents/Conduit Comsol/COMSOL/mc_failure.txt',' ',8);
r = data.data(:,1);
z = data.data(:,2);
[z,idx] = sort(z);
r = r(idx);


1 + 'hello'

plot(r,z,'DisplayName', 'COMSOL - no shear','LineWidth',3)

data = importdata('../Documents/Conduit Comsol/COMSOL/mc_failure_tau.txt',' ',8);
r = data.data(:,1);
z = data.data(:,2);
[z,idx] = sort(z);
r = r(idx);
plot(r,z,'DisplayName', 'COMSOL - with shear','LineWidth',3)
title('Failure Envelope Predicted by Elastic Analysis')
ylabel('z (m)')
xlabel('r (m)')


%% Match z to zvec
minz = min(z);
maxz = max(z);

zslice = zvec((zvec > minz) & (zvec < maxz));
pslice = pvec((zvec > minz) & (zvec < maxz));

rcomsol = nan(size(zslice));
for i = 1:length(rcomsol)
   [~,idx] = min(abs(z - zslice(i)));
   rcomsol(i) = r(idx);
end

%% with tau

%figure

bound_rt = @(p,s,tau) p*(3*s-p-1) - (s*(2*s-1) - 1/2*(tau)^2);
bound_rz = @(p,s,tau) (4*s-3*p-1)^2 - ((p-1)^2 + 4*tau^2);


s = 0:.01:1;
tau = 0.73;

prt = nan(size(s));
prz = nan(size(s));

for i = 1:length(s)
   prt(i) = fzero(@(x) bound_rt(x,s(i),tau*x),s(i)+tau);
   prz(i) = fzero(@(x) bound_rz(x,s(i),tau*x),2*s(i) - 1-tau);
end

plot(s,prt,'--r'); hold on
plot(s,prz,'--b');
ylim([0,1])
xlim([0,1])

%% failure surface
s = 0:0.01:1;

fail_rz = @(p,tau,c,q)   .5*(p+1+((p-1)^2 + 4*(tau*p)^2)^0.5) - c - q*0.5*(p+1-((p-1)^2 + 4*(tau*p)^2)^0.5);
fail_rt = @(p,s,tau,c,q) - 2*s + p + c + q*0.5*(p + 1 - ((p-1)^2 + 4*(tau*p)^2)^.5);
fail_tz = @(p,s,tau,c,q) 0.5*(p + 1 + ((p-1)^2 + 4*(tau*p)^2)^.5) - c - q*(2*s-p);

tau = 0.7;
c = 0;
phi = 38/180*pi;
%cohesion = 5e6; mu = tan(phi); C = 2*cohesion*((mu^2 + 1)^(1/2) + mu);

q = tan(pi/4 + 1/2*phi)^2;

frz = ones(size(s));
frz_lower = frz * fzero(@(x) fail_rz(x,tau,c,q), 0.237);
frz_upper = frz * fzero(@(x) fail_rz(x,tau,c,q), 1);
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

plot(s, frz_lower,'--r'); hold on;
plot(s, frz_upper,'--r');

plot(s, frt,'-b');
plot(s, ftz,'-g');
xlim([0,1])
ylim([0,1])
xlabel('S/\sigma_{zz}')
ylabel('p/\sigma_{zz}')

%%

c = 1;
l = nan(size(s));
for i = 1:length(s)
    l(i) = fail_rz(s(i),.7813,c,q1);
end
plot(s,l); hold on
plot(xlim,[0,0],'--r')
ylim([-2,2])

%% drucker prager analytical
syms sigz p lmda
A = diag([sigz, p, 2*lmda*sigz - p]);
t = trace(A)/3;
d = A - t;
simplify(sqrt(sum(d.*d,'all'))/sigz)
simplify(sqrt(1/6*(expand((sigz - p)^2 + (p-(2*lmda*sigz-p))^2 + (sigz-(2*lmda*sigz-p))^2))))
 


%% change of angle

pz = 0:.001:.99999;
th = nan(size(pz));
tp = 0.7813;
for i=1:length(pz)
    th(i) = 1/2*atand(2*tp*pz(i)/(pz(i)-1));
end

subplot(211)
plot(pz,th); hold on;
xlabel('p/\sigma_{zz}')
ylabel('rotation (deg)')

%% t-star

tstar = @(c,q)  ((1./((c + q).*(1-c))).^(1/2).*(c + q - 1))/2;

phi = 38*pi/180;
dt = 0.001;
c = 0:.01:1-dt;
q = tan(pi/4 + 1/2*phi)^2;

ts = tstar(c,q);
subplot(212)
semilogy(c,abs(ts)); hold on;
xlabel('c/\sigma_{zz}')
ylabel('\tau^*')

%% calculate the SRF
lambda = 0.5;
phi = 30*pi/180;

SRF = tan(phi)/(tan(2*(atan(sqrt(1/lambda))-pi/4)))

%% Aiy-Corona

SantaClara = [7, 9, 11, 14, 20, 24, 32, 37, 43, 45, 48, 66, 79, 91, 114, 138, 155, 175, 189, 196, 263, 302, 321, 375, 459, 542] ;