function [with_shear, no_shear, failure_shear, failure_no_shear] = plotfailureprofiles(A,Srr,Szz,Stt,Srz,zvec,pvec,plotfigs,porep)


tau = Srz./Szz;
pp = porep./Szz;
phi = A.mc.phi;
cohesion = A.mc.C(zvec); mu = tan(phi); c = 2*cohesion*((mu^2 + 1)^(1/2) + mu)./Szz;
C = 2*cohesion*((mu^2 + 1)^(1/2) + mu);

q = tan(pi/4 + 1/2*phi)^2;

fail_rz = @(p,tau,pp,c,q)   (.5*(p+1+((p-1)^2 + 4*(tau)^2)^0.5)-pp) - c - q*(0.5*(p+1-((p-1)^2 + 4*(tau)^2)^0.5)-pp);
fail_rt = @(p,s,tau,pp,c,q) -(2*s - p - pp) + c + q*(0.5*(p + 1 - ((p-1)^2 + 4*(tau)^2)^.5) - pp);
fail_tz = @(p,s,tau,pp,c,q) (0.5*(p + 1 + ((p-1)^2 + 4*(tau)^2)^.5) - pp) - c - q*(2*s-p - pp);

frz_low  = nan(length(zvec),1);
frz_high = nan(length(zvec),1);
frt = nan(length(zvec),1);
ftz = nan(length(zvec),1);

s = A.lambda; % Vertical gradient of horizontal stress

for i = 1:length(zvec)
    try
        frz_low(i) = fzero(@(x) fail_rz(x, tau(i), pp(i), c(i), q), 0.237);
    catch ME
    end
    try
        frz_high(i) = fzero(@(x) fail_rz(x, tau(i), pp(i), c(i), q), 5);
    catch ME
    end
    try
        frt(i) = fzero(@(x) fail_rt(x, s, tau(i), pp(i), c(i), q), [-1, 1]);
    catch ME
    end
    try
        ftz(i) = fzero(@(x) fail_tz(x, s, tau(i), pp(i), c(i), q), [-1, 1]);
    catch ME
    end
end

tstar = ((1./((c + q).*(1-c))).^(1/2).*(c + q - 1))/2;

% old code, not up to date with pore pressure


if plotfigs
lw =3; 
    
figure
clrs = parula(4);
subplot(121)
plot(Srz,zvec); hold on;
plot(Srr,zvec);
legend('\tau','p')
ylabel('z')
xlabel('Pa')

subplot(122)
plot(Srz./Srr,zvec);
hold on
plot(tstar, zvec, '--r')
ylabel('z')
xlabel('\tau/p')

figure

plotfrz = (frz_low>=frt);
plotfrt = ~plotfrz;
% add back in a point for overlap
plotfrz(min(find(plotfrt))) = 1;

p1 = plot(Srr./Szz, zvec,'-k','LineWidth',lw,'DisplayName','p/\sigma_{zz}'); hold on
xl = xlim;
p2 = plot(frz_low(plotfrz),zvec(plotfrz), 'Color', clrs(1,:),'LineWidth',lw, 'LineStyle', '-','DisplayName','r/z');
p3 = plot(frt(plotfrt),zvec(plotfrt),'Color',clrs(2,:),'LineWidth',lw, 'LineStyle', '-','DisplayName','r/\theta');
p4 = plot(ftz,zvec,'Color',clrs(3,:),'LineWidth',lw, 'LineStyle', '-','DisplayName','\theta/z');

% find C = szz
[~,ind] = min(abs(Szz-C+(q-1)*porep));

plot(xlim, [zvec(ind), zvec(ind)], '--r')


xlabel('p/\sigma_{zz}')
ylabel('z (m)');


xlim([0,1])
ylim([-6200,0])


x = [max(frz_low,frt);0];
y = [zvec;min(zvec)];
patch(x,y,'k','FaceAlpha',.1,'LineStyle', 'none');
idx = ~isnan(ftz);
x = [ftz(idx);1;1];
y = [zvec(idx);max(zvec(idx));min(zvec)];
patch(x,y,'k','FaceAlpha',.1,'LineStyle', 'none');

legend([p1, p2, p3, p4], 'Orientation','horizontal','Location','south')

end

C = 2*cohesion*((mu^2 + 1)^(1/2) + mu);

q = tan(pi/4 + 1/2*phi)^2;

rzrf = @(p, sigz, lambda, beta) sqrt((p./sigz - lambda)./(1/q - C./(q*sigz) - lambda));
rzrf2 = @(p,sigz, lambda, beta) (q*sigz.*(((lambda - p./sigz).*(C - sigz + lambda.*q.*sigz))./(q*sigz)).^(1/2))./(C - sigz + lambda.*q.*sigz);
rrtf = @(p, sig_z, x) (sig_z.*(q + 1).*(-((p./sig_z - x).*(C - sig_z.*x + q.*sig_z.*x))./(sig_z.*(q + 1))).^(1/2))./(C - sig_z.*x + q.*sig_z.*x);

rfailzr = rzrf(Srr, Szz, A.lambda);
rfailzr2 = rrtf(Srr,Szz,A.lambda);

if plotfigs
figure
plot(rfailzr*A.r, zvec, 'DisplayName', 'Analytical','LineWidth',3);
hold on;
plot(rfailzr2*A.r, zvec, 'DisplayName', 'Analytical2','LineWidth',3);

plot([A.r, A.r], ylim, '--r')

% xlim([100, 115]);
% ylim([-2900, -1100]);
% xlim([80, 180])
% ylim([-3000, 0])
legend show

end

[~,ind] = min(abs(zvec - A.fragdepth)); % find index of fragmentation depth
ind = ind - 1; % max shear stress is one step below

shear_condition = max(frz_low,frt);
no_shear_condition = max(rfailzr, rfailzr2);

with_shear = shear_condition(ind)/(Srr(ind)/Szz(ind)) - 1;

no_shear = no_shear_condition(ind) - 1;

failure_shear = any(shear_condition./(Srr./Szz) > 1);
failure_no_shear = any(no_shear_condition(ind) > 1);

end