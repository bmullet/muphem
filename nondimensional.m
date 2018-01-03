A = out{1}; zvec = out{2}; pvec = out{3}; ugvec = out{4}; umvec = out{5};
phivec = out{6}; rhogvec = out{7}; chidvec = out{8}; Qmvec = out{9}; Qgvec = out{10}; failure = out{11};


pc = max(pvec);
rhogc = pc/(A.Rw*A.T);
rhomc = A.rhom0;
uc = min(umvec);
delta = rhogc/rhomc;
L = A.depth;
deltapc = 1e6; % Note sure how to set this though


%Get non-dimensional variables
pnd = pvec/pc;
rhognd = rhogvec/rhogc;
rhomnd = ones(size(rhognd))*A.rhom0;
ugnd = ugvec/uc;
umnd = umvec/uc;
phi = phivec;
chid = chidvec;
znd = zvec/L;
rb = (3*phi/(4*pi*A.nb)).^(1/3);

%bubble growth phase
%% mass balance
Gp1 = phi + (1-phi).*(1/delta.*(rhomnd.*umnd./(rhognd.*ugnd)) - 1).*(A.hb.*chid./(1-chid));

c1 = Gp1/pnd;
c2 = (phi/ugnd);
c3 = (1-phi)/umnd;

subplot(1,3,1)
plot(c1,znd);
subplot(1,3,2)
plot(c2,znd);
subplot(1,3,3)
plot(c3,znd);

%% gas momentum
Ma = uc/sqrt(A.Rw*A.T);
Fr = uc/sqrt(A.g*A.r);
Stk = 2/9*rb.^2.*rhogc./(A.mu) ./ (A.r/uc);

c1 = phi/Ma^2 - (1-phi)*(1/delta).*(umnd.*ugnd.*rhomnd)./pnd.*(A.hb*chid./(1-chid));
c2 = phi.*rhognd.*ugnd;
c3 = -phi.*rhognd/Fr^2;
c4 = -phi.*(1-phi)./Stk.*(ugnd - umnd);
subplot(1,4,1)
plot(-1*c1,znd);
set(gca,'XScale','log')
ylim([0 1])
subplot(1,4,2)
plot(c2,znd);
subplot(1,4,3)
plot(c3,znd);
subplot(1,4,4)
plot(abs(c4),znd);
set(gca,'XScale','log')
ylim([0 1])
xlim([0 100])

figure
subplot(1,2,1)
plot(phi.*(1-phi)./Stk,znd)
subplot(1,2,2)
plot(ugnd - umnd, znd)

%% melt momentum
Ma = uc/sqrt(A.Rw*A.T);
Fr = uc/sqrt(A.g*A.r);
Stk = 2/9*rb.^2.*rhogc./(A.mu) ./ (A.r/uc);
Re = A.r*rhomc*uc/(8*A.mu);

c1 = (1-phi)*delta/Ma^2 + (1-phi).*(umnd.*umnd.*rhomnd)./pnd.*(A.hb*chid./(1-chid));
c2 = (1-phi).*rhomnd.*umnd;
c3 = -(1-phi).*rhomnd/Fr^2;
c4 = phi.*(1-phi)./Stk.*(ugnd - umnd);
c5 = umnd/Re;

subplot(1,5,1)
plot(-1*c1,znd);
set(gca,'XScale','log')
ylim([0 1])
subplot(1,5,2)
plot(c2,znd);
subplot(1,5,3)
plot(c3,znd);
subplot(1,5,4)
plot(c4,znd);
ylim([0 1])
subplot(1,5,5)
plot(c5,znd)