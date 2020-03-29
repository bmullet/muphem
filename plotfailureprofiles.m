function [with_shear, no_shear, failure] = plotfailureprofiles(A,Srr,Szz,Stt,Srz,zvec,pvec,plotfigs)

close all;

tau = Srz./Szz;
phi = A.mc.phi;
cohesion = A.mc.C; mu = tan(phi); c = 2*cohesion*((mu^2 + 1)^(1/2) + mu)./Szz;

q = tan(pi/4 + 1/2*phi)^2;

fail_rz = @(p,tau,c,q)   .5*(p+1+((p-1)^2 + 4*(tau)^2)^0.5) - c - q*0.5*(p+1-((p-1)^2 + 4*(tau)^2)^0.5);
fail_rt = @(p,s,tau,c,q) -2*s + p + c + q*0.5*(p + 1 - ((p-1)^2 + 4*(tau)^2)^.5);
fail_tz = @(p,s,tau,c,q) 0.5*(p + 1 + ((p-1)^2 + 4*(tau)^2)^.5) - c - q*(2*s-p);

frz_low  = nan(length(zvec),1);
frz_high = nan(length(zvec),1);
frt = nan(length(zvec),1);
ftz = nan(length(zvec),1);

s = A.lambda; % Vertical gradient of horizontal stress

for i = 1:length(zvec)
    try
        frz_low(i) = fzero(@(x) fail_rz(x, tau(i), c(i), q), 0.237);
    catch ME
    end
    try
        frz_high(i) = fzero(@(x) fail_rz(x, tau(i), c(i), q), 5);
    catch ME
    end
    try
        frt(i) = fzero(@(x) fail_rt(x, s, tau(i), c(i), q), [-1, 1]);
    catch ME
    end
    try
        ftz(i) = fzero(@(x) fail_tz(x, s, tau(i), c(i), q), [-1, 1]);
    catch ME
    end
end

tstar = ((1./((c + q).*(1-c))).^(1/2).*(c + q - 1))/2;

tfail = ((c + pvec./Szz*q - 1).*(c - pvec./Szz + q)).^(1/2)./(pvec./Szz + pvec./Szz*q);

if plotfigs
figure
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
plot(Srr./Szz, zvec,'DisplayName','p/\sigma_{zz})'); hold on
xl = xlim;
plot(frz_low,zvec,'DisplayName','r/z');
plot(frt,zvec,'DisplayName','r/\theta');
plot(ftz,zvec,'DisplayName','\theta/z');
xlabel('p/S_{zz}')
ylabel('z');
xlim(xl)
legend()


figure
end

C = 2*cohesion*((mu^2 + 1)^(1/2) + mu);

q = tan(pi/4 + 1/2*phi)^2;

rzrf = @(p, sigz, lambda, beta) sqrt((p./sigz - lambda)./(1/q - C./(q*sigz) - lambda));

rfailzr = rzrf(Srr, Szz, A.lambda);

if plotfigs
plot(rfailzr*A.r, zvec, 'DisplayName', 'Analytical','LineWidth',3);
hold on;

plot([A.r, A.r], ylim, '--r')

% xlim([100, 115]);
% ylim([-2900, -1100]);
% xlim([80, 180])
% ylim([-3000, 0])
legend show

end

[~,ind] = min(abs(zvec - A.fragdepth)); % find index of fragmentation depth
ind = ind - 1; % max shear stress is one step below

with_shear = frz_low(ind)/(Srr(ind)/Szz(ind)) - 1;
no_shear = rfailzr(ind) - 1;

failure = any(frz_low./(Srr./Szz) > 1);


end