function [ ] = plotmuphem(A,z,p,ug,um,phi,rhog,chid)
% This code is designed to plot the results of the multiphase code.

set(0, 'DefaultAxesFontWeight', 'normal', ...
      'DefaultAxesFontSize', 20, ...
      'DefaultAxesTitleFontWeight', 'bold', ...
      'DefaultAxesTitleFontSizeMultiplier', 1, ...
      'DefaultLineLineWidth',3) ;
% This is a new comment
close all
% A.delF = 0;
% eos = eosf(0);
% 
% try
%     [rhog, chid, phi] = eos.calcvars(A,um,p);
%     A.delF = 1;
%     eos = eosf(1);
%     [rhog2, chid2, phi2] = eos.calcvars(A,um,p);
%     rhog(z<A.fragdepth) = rhog2(z<A.fragdepth);
%     chid(z<A.fragdepth) = chid2(z<A.fragdepth);
%     phi(z<A.fragdepth) = phi2(z<A.fragdepth);
%     phi(z<A.exdepth) = 0;
% catch
%     disp('Note: This solution does not fragment')
%     A.delF = 1;
%     eos = eosf(1);
%     [rhog, chid, phi] = eos.calcvars(A,um,p);
%     phi(z<A.exdepth) = 0;
% end

rhom = A.rhom0;
Qg = phi.*rhog.*ug;
Qm = (1-phi).*rhom.*um;
Qt = Qg+Qm;
pplot = p/1e6;
[~,i] = min(abs(p-A.Pcrit*0.99));
z = abs(z);
zprint = z/1e3;
zex = z(i)*ones(size(p))/1e3;

disp('Pressure at conduit exit:')
disp(min(p));

figure()
% um, ug vs. z
hold on
subplot(1,8,1)
plot(um,zprint,ug,zprint)
hold on
set(gca,'XScale','log')
xlabel('Velocity (m/s)')
set(gca,'Ydir','reverse')
ylabel('depth (km)')


% Pressure v. z
subplot(1,8,2)
hold on
plot(pplot,zprint)
xlabel('Pressure (MPa)')
set(gca,'Ydir','reverse')
ylabel('depth (km)')

% Phi v. z
subplot(1,8,3)
hold on
plot(phi,zprint)
xlabel('\phi')
set(gca,'Ydir','reverse')
ylabel('depth (km)')

% Q
subplot(1,8,4)
hold on
plot(Qg./Qm,zprint)
xlabel('Q_g/Q_m (kg/m^2-s)')
set(gca,'Ydir','reverse')
ylabel('depth (km)')

% chi v. z
subplot(1,8,5)
hold on
plot(chid,zprint)
xlabel('\chi_d')
set(gca,'Ydir','reverse')
ylabel('depth (km) ')

% rhog v. z
subplot(1,8,6)
hold on
plot(rhog,zprint)
xlabel('\rho_g (kg/m^3)')
set(gca,'Ydir','reverse')
ylabel('depth (km) ')

% mu v. z
subplot(187)
hold on
mu = A.mu(phi,p);
plot(log10(mu(zprint>abs(A.fragdepth)/1e3)),zprint(zprint>abs(A.fragdepth/1e3)))
xlabel('log_{10} \mu')
set(gca,'Ydir','reverse')
ylabel('depth (km)')

% chi_c v. z
subplot(188)
hold on
xc = A.xc(p);
plot(xc, zprint)
xlabel('\chi_c')
set(gca,'Ydir','reverse')
ylabel('depth (km)')

for i=1:8
  subplot(1,8,i)
  xl=xlim;
  hold on
  plot([min(xl) max(xl)],[zex zex],'k--')
  plot([min(xl) max(xl)],[abs(A.fragdepth)/1e3, abs(A.fragdepth)/1e3],'r--')
  box('on')
  hold off
end

set(gcf,'Units','inches',...
 'Position',[0 0 20 10])

subplot(1,8,1)
legend('Melt', 'Gas')
legend('show')
drawnow

% % Solubility Law
% figure
% hold on
% subplot(1,2,1)
% title('Solubility v. z')
% hold on
% plot(chid,z)
% plot(chid,zex,'--r')
% xlabel('\chi_d')
% ylabel('z (m)')
% subplot(1,2,2)
% title('Solubility v. P')
% hold on
% plot(p,chid)
% xlabel('p')
% ylabel('\chi_d')
% 
% u vs. P

figure
%subplot(1,2,1)
plot(p,um,p,ug)
xlabel('P (Pa)')
ylabel('u (m/s)')
legend('Melt','Gas')
set(gca,'XScale','log')
set(gca,'YScale','log')

figure
subplot(1,3,1)
plot(Qm+Qg,z)
set(gca,'Ydir','reverse')
hold on
xlabel('Qt')
ylabel('z (m)')
subplot(1,3,2)
plot(Qm,z)
set(gca,'Ydir','reverse')
xlabel('Qm')
ylabel('z (m)')
subplot(1,3,3)
plot(Qg,z)
set(gca,'Ydir','reverse')
xlabel('Qg')
ylabel('z (m)')

% figure
% subplot(1,3,1)
% plot(um,z)
% xlabel('um')
% subplot(1,3,2)
% plot(phi,z)
% xlabel('phi')
% subplot(1,3,3)
% plot(rhom,z)
% xlabel('rhom')
% 
% figure
% plot(phi,z)
% title('phi')
% 
% figure
% plot((1-phi(z<A.fragdepth)).*um(z<A.fragdepth),z(z<A.fragdepth))
% title('Qm')
% 

ed = A.exdepth;
fd = A.fragdepth;
bgi = (z>=ed&z<=fd);
A.delF = 1;
delF = 1;
RHS1 = zeros(nnz(bgi),3);
phibg = phi(bgi);
umbg = um(bgi);
ugbg = ug(bgi);
rhogbg = rhog(bgi);
pbg = p(bgi);
Qmbg = Qm(bgi); chidbg = chid(bgi);
Fmg = zeros(nnz(bgi),1);
F1 = Fmg;
rbvec = Fmg;
Is = Fmg;
Fw = Fmg;
% for i=1:nnz(bgi)
% rb = (3*phibg(i)/(4*pi*A.nb))^(1/3); %bubble radius
% rbvec(i) = 1/rb^2;
% Fmg(i) = -9/2*phibg(i).*(1-phibg(i)).*A.mu.*(ugbg(i)-umbg(i))/rb^2;
% F1(i) = -phibg(i)*rhogbg(i)*A.g;
% Is(i) = -Qmbg(i)*A.hs*A.hb*pbg(i)^(A.hb-1)/(1-chidbg(i));
% RHS1(i,:) = twophaseODE(z(i),[p(i), ug(i), um(i)],A);
% Fw(i) = 8*A.mu*umbg(i)/A.r^2;
% end

figure()
plot(ugbg./umbg,z(bgi))
xlabel('u_g/u_m')
ylabel('z (m)')

% 
% figure()
% subplot(2,2,1)
% plot(Is,z(bgi)); title('I^*')
% subplot(2,2,2)
% plot(F1,z(bgi)); title('-\phi \rho_g g')
% subplot(2,2,3)
% plot(Fmg,z(bgi)); title('Interphase')
% subplot(2,2,4)
% plot(Fw,z(bgi)); title('Wall Friction')
% 
% ed = A.exdepth;
% fd = A.fragdepth;
% fi = (z>=fd);
% A.delF = 1;
% delF = 1;
% phif = phi(fi);
% umf = um(fi);
% ugf = ug(fi);
% rhogf = rhog(fi);
% pf = p(fi);
% Qmf = Qm(fi); chidbg = chid(fi);
% Fmg = zeros(nnz(fi),1);
% F1 = Fmg;
% Is = Fmg;
% Fw = Fmg;
% for i=1:nnz(fi)
% Fmg(i) = -3/8*phif(i)*(1-phif(i))*A.dragC/A.Rash*rhogf(i)*abs(ugf(i)-umf(i))*(ugf(i)-umf(i));
% F1(i) = -phif(i)*rhogf(i)*A.g;
% %RHS1(i,:) = twophaseODE(z(i),[p(i), ug(i), um(i)],A);
% tau = A.f0.*0.5*rhogf(i).*ugf(i).^2;
% Fw(i) = 2*tau./A.r;
% end
% A.delF = 0;

% figure()
% subplot(2,2,2)
% plot(F1,z(fi)); title('-\phi \rho_g g')
% subplot(2,2,3)
% plot(Fmg,z(fi)); title('Interphase')
% subplot(2,2,4)
% plot(Fw,z(fi)); title('Wall Friction')
% 

%RHS2 = twophaseODE(z(z>=fd),[p(z>=fd), ug(z>=fd), um(z>=fd)],A);
%RHS = [RHS1;RHS2];

%plot(RHS,z);

% subplot(1,2,2)
% plot(pplot,phi)
% xlabel('P (MPa)')
% ylabel('\phi')
% set(gca,'XScale','log')
% set(gca,'YScale','log')
ppore = A.Patm_ + 10000*(A.depth-z);
ppore = A.rhow*A.g*(A.depth-z)+A.Patm_;
plith = A.k.rho*A.g*(A.depth-z)+A.Patm_;
figure
plot(ppore,z,p,z,plith,z)
xlabel('Pressure (Pa)')
ylabel('Depth (z)')
legend('Pore Pressure','Magma Pressure','Lithostatic')

figure
plot(um,z,ug,z)
hold on
set(gca,'XScale','log')
xlabel('Velocity (m/s)')
ylabel('z (m)')
xl=xlim;

plot([min(xl) max(xl)],[-A.fragdepth -A.fragdepth],'r--')
set(gca,'Ydir','reverse')

figure
semilogx(p,z);
set(gca,'XScale','log')
xlabel('Pressure (Pa)')
set(gca,'Ydir','reverse')
ylabel('depth (m)')

figure
loglog(p,phi);
set(gca,'XScale','log')
xlabel('Pressure (Pa)')
ylabel('\phi')


% plot elongation strain rate
k = 0.01;
Ginf = 20e9;
ez = DGradient(um,z);
figure
plot(-ez(z>-A.fragdepth),z(z>-A.fragdepth));
hold on;
crit = @(Ginf) k*Ginf./mu(z>-A.fragdepth);
plot(crit(3e9), z(z>-A.fragdepth));
plot(crit(15e9), z(z>-A.fragdepth));
plot(crit(30e9), z(z>-A.fragdepth));
legend('Strain rate','Critical - 3GPa','Critical - 15GPa','Critical - 30GPa')
xlabel('\epsilon_{zz}')
set(gca,'Ydir','reverse')
ylabel('z')



