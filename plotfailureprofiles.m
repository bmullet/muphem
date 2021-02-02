function [with_shear_low, no_shear_low, with_shear_high, no_shear_high, failure_shear, failure_no_shear, failure_mechanisms] = plotfailureprofiles(A,Srr,Szz,~,Srz,zvec,~,plotfigs,porep)


tau = Srz./Szz;
taup = Srz./Srr; % normalized by pressure
pp = porep./Szz;
phi = A.mc.phi;
cohesion = A.mc.C(zvec); mu = tan(phi); c = 2*cohesion*((mu^2 + 1)^(1/2) + mu)./Szz;
C = 2*cohesion*((mu^2 + 1)^(1/2) + mu);

q = tan(pi/4 + 1/2*phi)^2;

% These expressions are derived analytically in the derivations.m script
% TODO: the expressions have additional components that might come into
% play depending on taup. Add those.
fail_rz_low = @(taup,pp,c,q)   -(c + pp - c*q - 2*pp*q + pp*q^2 + q*(4*c^2*taup^2 + c^2 - 8*c*pp*q*taup^2 - 2*c*pp*q + 8*c*pp*taup^2 + 2*c*pp + 4*c*q*taup^2 + 2*c*q - 4*c*taup^2 - 2*c + 4*pp^2*q^2*taup^2 + pp^2*q^2 - 8*pp^2*q*taup^2 - 2*pp^2*q + 4*pp^2*taup^2 + pp^2 - 4*pp*q^2*taup^2 - 2*pp*q^2 + 8*pp*q*taup^2 + 4*pp*q - 4*pp*taup^2 - 2*pp + q^2 - 4*q*taup^2 - 2*q + 1)^(1/2) - q^2 + (4*c^2*taup^2 + c^2 - 8*c*pp*q*taup^2 - 2*c*pp*q + 8*c*pp*taup^2 + 2*c*pp + 4*c*q*taup^2 + 2*c*q - 4*c*taup^2 - 2*c + 4*pp^2*q^2*taup^2 + pp^2*q^2 - 8*pp^2*q*taup^2 - 2*pp^2*q + 4*pp^2*taup^2 + pp^2 - 4*pp*q^2*taup^2 - 2*pp*q^2 + 8*pp*q*taup^2 + 4*pp*q - 4*pp*taup^2 - 2*pp + q^2 - 4*q*taup^2 - 2*q + 1)^(1/2) - 1)/(2*(q^2*taup^2 + 2*q*taup^2 + q + taup^2));
fail_rz_high = @(taup,pp,c,q) (c*q - pp - c + 2*pp*q - pp*q^2 + q*(4*c^2*taup^2 + c^2 - 8*c*pp*q*taup^2 - 2*c*pp*q + 8*c*pp*taup^2 + 2*c*pp + 4*c*q*taup^2 + 2*c*q - 4*c*taup^2 - 2*c + 4*pp^2*q^2*taup^2 + pp^2*q^2 - 8*pp^2*q*taup^2 - 2*pp^2*q + 4*pp^2*taup^2 + pp^2 - 4*pp*q^2*taup^2 - 2*pp*q^2 + 8*pp*q*taup^2 + 4*pp*q - 4*pp*taup^2 - 2*pp + q^2 - 4*q*taup^2 - 2*q + 1).^(1/2) + q^2 + (4*c^2*taup^2 + c^2 - 8*c*pp*q*taup^2 - 2*c*pp*q + 8*c*pp*taup^2 + 2*c*pp + 4*c*q*taup^2 + 2*c*q - 4*c*taup^2 - 2*c + 4*pp^2*q^2*taup^2 + pp^2*q^2 - 8*pp^2*q*taup^2 - 2*pp^2*q + 4*pp^2*taup^2 + pp^2 - 4*pp*q^2*taup^2 - 2*pp*q^2 + 8*pp*q*taup^2 + 4*pp*q - 4*pp*taup^2 - 2*pp + q^2 - 4*q*taup^2 - 2*q + 1).^(1/2) + 1)/(2*(q^2*taup^2 + 2*q*taup^2 + q + taup^2));
fail_rt_right = @(k,taup,pp,c,q) -(2*c - 4*k + 2*pp + q + c*q - 2*k*q - pp*q - q*(4*c^2*taup^2 + c^2 - 16*c*k*taup^2 - 4*c*k - 8*c*pp*q*taup^2 - 2*c*pp*q + 8*c*pp*taup^2 + 2*c*pp + 4*c*q*taup^2 + 2*c*q + 2*c + 16*k.^2*taup^2 + 4*k.^2 + 16*k*pp*q*taup^2 + 4*k*pp*q - 16*k*pp*taup^2 - 4*k*pp - 8*k*q*taup^2 - 4*k*q - 4*k + 4*pp^2*q^2*taup^2 + pp^2*q^2 - 8*pp^2*q*taup^2 - 2*pp^2*q + 4*pp^2*taup^2 + pp^2 - 4*pp*q^2*taup^2 - 2*pp*q^2 + 4*pp*q*taup^2 + 2*pp + q^2 + 2*q + 1).^(1/2) - pp*q^2 + q^2)./(2*(- q^2*taup^2 + q + 1));
fail_tz_left = @(k,taup,pp,c,q) (c + pp - q + 2*c*q + 2*k*q + pp*q - (4*c^2*taup^2 + c^2 + 16*c*k*q*taup^2 + 4*c*k*q - 8*c*pp*q*taup^2 - 2*c*pp*q + 8*c*pp*taup^2 + 2*c*pp - 2*c*q - 4*c*taup^2 - 2*c + 16*k.^2*q^2*taup^2 + 4*k.^2*q^2 - 16*k*pp*q^2*taup^2 - 4*k*pp*q^2 + 16*k*pp*q*taup^2 + 4*k*pp*q - 4*k*q^2 - 8*k*q*taup^2 - 4*k*q + 4*pp^2*q^2*taup^2 + pp^2*q^2 - 8*pp^2*q*taup^2 - 2*pp^2*q + 4*pp^2*taup^2 + pp^2 + 2*pp*q^2 + 4*pp*q*taup^2 - 4*pp*taup^2 - 2*pp + q^2 + 2*q + 1).^(1/2) + 4*k*q^2 - 2*pp*q^2 - 1)./(2*(q^2 + q - taup^2));

frz_low  = nan(length(zvec),1);
frz_high = nan(length(zvec),1);
frt_right = nan(length(zvec),1);
ftz_left = nan(length(zvec),1);

frz_low_notau  = nan(length(zvec),1);
frz_high_notau = nan(length(zvec),1);
frt_right_notau = nan(length(zvec),1);
ftz_left_notau = nan(length(zvec),1);

% TODO: change these to analytical expressions
s = A.lambda(zvec); % Vertical gradient of horizontal stress

for i = 1:length(zvec)
    try
        frz_low(i) = fail_rz_low(taup(i), pp(i), c(i), q);
    catch ME
    end
    try
        frz_high(i) = fail_rz_high(taup(i), pp(i), c(i), q);
    catch ME
    end
    try
        frt_right(i) = fail_rt_right(s(i), tau(i), pp(i), c(i), q);
    catch ME
    end
    try
        ftz_left(i) = fail_tz_left(s(i), tau(i), pp(i), c(i), q);
    catch ME
    end
    try
        frz_low_notau(i) = fail_rz_low(0, pp(i), c(i), q);
    catch ME
    end
    try
        frz_high_notau(i) = fail_rz_high(0, pp(i), c(i), q);
    catch ME
    end
    try
        frt_right_notau(i) = fail_rt_right(s(i), 0, pp(i), c(i), q);
    catch ME
    end
    try
        ftz_left_notau(i) = fail_tz_left(s(i), 0, pp(i), c(i), q);
    catch ME
    end
end

tstar = ((1./((c + q).*(1-c))).^(1/2).*(c + q - 1))/2;



if plotfigs
    lw =3;
    
    % Figure for tau and p
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
    
    % Failure figure
    figure
     
    plotfrz = (frz_low>=frt_right);
    plotfrt = ~plotfrz;
    plotfrt(min(find(plotfrt)):end) = 1;
    plotfrz = ~plotfrt;
    % add back in a point for overlap
    plotfrz(min(find(plotfrt))) = 1;
    
    plotftz = (frz_high >= ftz_left);
    plotfrz_high = ~plotftz;
    plotfrz_high(min(find(plotftz))) = 1;
    
    p1 = plot(Srr./Szz, zvec,'-k','LineWidth',lw,'DisplayName','p/\sigma_{zz}'); hold on
    xl = xlim;
    p2 = plot(frz_low(plotfrz),zvec(plotfrz), 'Color', clrs(1,:),'LineWidth',lw, 'LineStyle', '-','DisplayName','r/z');
    p3 = plot(frt_right(plotfrt),zvec(plotfrt),'Color',clrs(2,:),'LineWidth',lw, 'LineStyle', '-','DisplayName','r/\theta');
    p4 = plot(ftz_left(plotftz),zvec(plotftz),'Color',clrs(3,:),'LineWidth',lw, 'LineStyle', '-','DisplayName','\theta/z');
    p5 = plot(frz_high(plotfrz_high),zvec(plotfrz_high),'Color',clrs(4,:),'LineWidth',lw, 'LineStyle', '-','DisplayName','z/r');
    
    % find C = szz
    %[~,ind] = min(abs(Szz-C+(q-1)*porep));
    
    %plot(xlim, [zvec(ind), zvec(ind)], '--r')
    
    plot(xlim, [A.fragdepth, A.fragdepth], '--r','LineWidth',1)
    
    xlabel('p/\sigma_{zz}')
    ylabel('z (m)');
    
    
    xlim([0,4])
    ylim([-6200,0])
    
    
    x = [max(frz_low,frt_right);0];
    y = [zvec;min(zvec)];
    patch(x,y,'k','FaceAlpha',.1,'LineStyle', 'none');
    idx = ~isnan(ftz_left);
    x = [ftz_left(idx);1;1];
    y = [zvec(idx);max(zvec(idx));min(zvec)];
    patch(x,y,'k','FaceAlpha',.1,'LineStyle', 'none');
    
    legend([p1, p2, p3, p4, p5], 'Orientation','horizontal','Location','south')
    
    % Frictional faulting profile
    
    
end

% Get failure conditions

[with_shear_low, with_shear_low_mech] = min([min(real(Srr./Szz - frz_low)), min(real(Srr./Szz - frt_right))]);
[with_shear_high, with_shear_high_mech] = min([min(real(frz_high - Srr./Szz)), min(real(ftz_left - Srr./Szz))]);
[no_shear_low, no_shear_low_mech] = min([min(real(Srr./Szz - frz_low_notau)), min(real(Srr./Szz - frt_right_notau))]);
[no_shear_high, no_shear_high_mech] = min([min(real(frz_high_notau - Srr./Szz)), min(real(ftz_left_notau - Srr./Szz))]);


% Keep which failure mechanism is responsible
pz = Srr./Szz;
if with_shear_low_mech == 0
    with_shear_low_mech = "r/z"; 
else
    [~,I] = min(Srr./Szz - frt_right);
    if pz(I) < 1
        with_shear_low_mech = "r/t";
    else
        with_shear_low_mech = "z/t";
    end
end

if with_shear_high_mech == 0
    [~,I] = min(frz_high - Srr./Szz);
    if pz(I) < 1    
        with_shear_high_mech = "r/z";
    else
        with_shear_high_mech = "z/r";
    end
else
    [~,I] = min(ftz_left - Srr./Szz);
    if pz(I) < 1
        with_shear_high_mech = "t/z";
    else
        with_shear_high_mech = "t/r";
    end
end

if no_shear_low_mech == 0
    no_shear_low_mech = "r/z";
else
    [~,I] = min(Srr./Szz - frt_right_notau);
    if pz(I) < 1
        no_shear_low_mech = "r/t";
    else
        no_shear_low_mech = "z/t";
    end
end

if no_shear_high_mech == 0
    [~,I] = min(frz_high_notau - Srr./Szz);
    if pz(I) < 1    
        no_shear_high_mech = "r/z";
    else
        no_shear_high_mech = "z/r";
    end
else
    [~,I] = min(ftz_left_notau - Srr./Szz);
    if pz(I) < 1
        no_shear_high_mech = "t/z";
    else
        no_shear_high_mech = "t/r";
    end
end

failure_shear = (with_shear_low < 0) || (with_shear_high < 0);
failure_no_shear = (no_shear_low < 0) || (no_shear_high < 0);

failure_mechanisms = {with_shear_low_mech, with_shear_high_mech, no_shear_low_mech, no_shear_high_mech};
end