function [with_shear, no_shear, failure_shear, failure_no_shear] = plotfailureprofiles(A,Srr,Szz,~,Srz,zvec,~,plotfigs,porep)


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

frz_low_notau  = nan(length(zvec),1);
frz_high_notau = nan(length(zvec),1);
frt_notau = nan(length(zvec),1);
ftz_notau = nan(length(zvec),1);

% TODO: change these to analytical expressions
s = A.lambda(zvec); % Vertical gradient of horizontal stress

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
        frt(i) = fzero(@(x) fail_rt(x, s(i), tau(i), pp(i), c(i), q), [-1, 1]);
    catch ME
    end
    try
        ftz(i) = fzero(@(x) fail_tz(x, s(i), tau(i), pp(i), c(i), q), [-1, 1]);
    catch ME
    end
    try
        frz_low_notau(i) = fzero(@(x) fail_rz(x, 0, pp(i), c(i), q), 0.237);
    catch ME
    end
    try
        frz_high_notau(i) = fzero(@(x) fail_rz(x, 0, pp(i), c(i), q), 5);
    catch ME
    end
    try
        frt_notau(i) = fzero(@(x) fail_rt(x, s(i), 0, pp(i), c(i), q), [-1, 1]);
    catch ME
    end
    try
        ftz_notau(i) = fzero(@(x) fail_tz(x, s(i), 0, pp(i), c(i), q), [-1, 1]);
    catch ME
    end
end

tstar = ((1./((c + q).*(1-c))).^(1/2).*(c + q - 1))/2;



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
    plotfrt(min(find(plotfrt)):end) = 1;
    plotfrz = ~plotfrt;
    % add back in a point for overlap
    plotfrz(min(find(plotfrt))) = 1;
    
    p1 = plot(Srr./Szz, zvec,'-k','LineWidth',lw,'DisplayName','p/\sigma_{zz}'); hold on
    xl = xlim;
    p2 = plot(frz_low(plotfrz),zvec(plotfrz), 'Color', clrs(1,:),'LineWidth',lw, 'LineStyle', '-','DisplayName','r/z');
    p3 = plot(frt(plotfrt),zvec(plotfrt),'Color',clrs(2,:),'LineWidth',lw, 'LineStyle', '-','DisplayName','r/\theta');
    p4 = plot(ftz,zvec,'Color',clrs(3,:),'LineWidth',lw, 'LineStyle', '-','DisplayName','\theta/z');
    
    % find C = szz
    %[~,ind] = min(abs(Szz-C+(q-1)*porep));
    
    %plot(xlim, [zvec(ind), zvec(ind)], '--r')
    
    plot(xlim, [A.fragdepth, A.fragdepth], '--r','LineWidth',1)
    
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
    
    % Frictional faulting profile
    
    
end


[~,ind] = min(abs(zvec - A.fragdepth)); % find index of fragmentation depth
ind = ind - 1; % max shear stress is one step below

shear_condition = max(frz_low,frt);
no_shear_condition = max(frz_low_notau, frt_notau);

with_shear = shear_condition(ind)/(Srr(ind)/Szz(ind)) - 1;
no_shear = no_shear_condition(ind)/(Srr(ind)/Szz(ind)) - 1;

failure_shear = any(shear_condition./(Srr./Szz) > 1);
failure_no_shear = any(no_shear_condition./(Srr./Szz) > 1);

end