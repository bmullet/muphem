function [] = plotbalances(A,z,p,ug,um,phi,rhog,chid)
%%% ASSUMES NO CRYSTAL GROWTH %%%


% Get indecies for pre and post fragmentation
prefrag = (z<=A.fragdepth);
postfrag = (z>A.fragdepth);

% Get derivatives and other values for balance equations
Qm = um.*(1-phi).*A.rhom0;
Qg = ug.*phi.*rhog;
du = ug-um;
dum = DGradient(um,z);
dug = DGradient(ug,z);
dphi = DGradient(phi,z);
drhog = DGradient(rhog,z);
ddu = DGradient((ug-um),z);
dp = DGradient(p,z);
rhom = A.rhom0;

[I, beta] = calcI(A, Qm, p, dp, prefrag, postfrag);
Fmg = interphase((ug-um), p, phi, rhog, A)/5.7e11;

figure
%Gas momentum balance
inertia = phi.*rhog.*ug.*dug;
pgrad = phi.*dp;
gravity = -phi.*rhog.*9.8;
mass = I.*ug;
Fgw = calcFgw(A,rhog,ug,prefrag);


plot(inertia,z); hold on;
plot(pgrad,z);
plot(gravity,z);
plot(-Fmg,z);
plot(mass,z);
plot(-Fgw,z);
xlim([-1e4,1e4]);
title('Gas momentum')

plot(inertia  + pgrad - gravity - Fmg - mass + Fgw,z,'--r')
legend('Inertia','Pressure','Gravity','Interphase','Mass','Wall','Balance');

%Melt momentum balance
% figure
% hold on;
% inertia = (1-phi).*rhom.*um.*dum;
% pgrad = (1-phi).*dp;
% gravity = -(1-phi).*rhom.*9.8;
% mass = I.*um;
% Fmw = calcFmw(A,um,phi,p,postfrag);
% 
% title('Melt momentum')
% plot(inertia,z); hold on;
% plot(pgrad,z);
% plot(-gravity,z);
% plot(Fmg,z);
% plot(mass,z);
% plot(Fmw,z);
% 
% 
% plot(inertia + pgrad  - gravity + Fmg - mass + Fmw,z,'--r')
% 
% legend('Inertia','Pressure','Gravity','Interphase','Mass','Wall','Balance');
% 
% 
% figure
% % Gas mass
% title('Gas mass')
% plot(phi.*rhog.*dug,z); hold on;
% plot(dphi.*rhog.*ug,z);
% plot(phi.*drhog.*ug,z);
% plot(-I,z);
% plot(phi.*rhog.*dug + dphi.*rhog.*ug + phi.*drhog.*ug - I,z, '--r');
% legend('dug','dphi','drhog','I', 'balance')
% xlim([-.1,.1])
% 
% figure
% % Melt mass
% title('Melt mass')
% plot((1-phi).*rhom.*dum,z); hold on;
% plot(-dphi.*rhom.*um,z);
% plot(I,z);
% plot((1-phi).*rhom.*dum -dphi.*rhom.*um + I,z, '--r');
% legend('dum','dphi','I', 'balance')
% xlim([-6,6])
% 
% 
% % Momentum 1
% % figure
% % alpha = -(Qm+Qg)./(p.*(1./ug + (1-phi)./phi.*1./um));
% % gamma = Qm./Qg.*(rhog.*ug./(rhom.*um) - 1).*beta;
% % gammat = 1 + alpha.*(1-gamma) + Qm.*beta.*du;
% % md = -(alpha.*p./ug + Qg);
% % 
% % 
% % pgrad = gammat.*dp;
% % deltau = -md.*ddu;
% % gravity = -(phi.*rhog + (1-phi).*rhom)*9.8;
% % plot(pgrad,z); hold on;
% % plot(deltau,z);
% % plot(gravity,z);
% % plot(Fgw,z);
% % plot(Fmw,z);
% % plot(pgrad - gravity + Fmw + Fgw,z,'--r');
% % 
% % legend('Pressure','\Delta u','Gravity','Fgw','Fmw','balance')
% % 




end

function [I,beta] = calcI(A,Qm,p,dp,~,postfrag)
    eos = eosf(1); % pre-fragmentation
    
    beta = eos.calcbeta(A,p./A.Pchamber);
    I = -Qm.*beta.*dp/A.Pchamber;
    
    I(postfrag) = 0;
    beta(postfrag) = 0;
    
    
end

function [Fgw] = calcFgw(A,rhog,ug,prefrag)

Fgw = A.f0.*rhog.*ug.^2./A.r;
Fgw(prefrag) = 0;

end

function [Fmw] = calcFmw(A,um,phi,p,postfrag)

Fmw = 8*A.mu(phi,p).*um./A.r^2;
Fmw(postfrag) = 0;

end


