function [] = plot_failure(z,pvec,phivec,A,Srr,Szz,Stt,Smax,Sfail,op)
%     Smin = min([Srr Szz Stt],[],2);
%     Smax = max([Srr Szz Stt],[],2); 
%     
%     qu = 2*A.mc.C*tan(deg2rad(45)+A.mc.phi/2);
%     Sfail = qu + Smin*tan(deg2rad(45)+A.mc.phi/2)^2;
%     
    % Get fragmentation depth
    z = (A.depth-z)/1e3;
    Zfr = z(find(phivec>A.phi0,1,'first'));
    Srrplot = Srr/1e6;
    Szzplot = Szz/1e6;
    Sttplot = Stt/1e6;
    pplot = pvec/1e6;
    figure()

    subplot(1,2,1);
    plot(Srrplot,z,Szzplot,z,Sttplot,z,pplot,z);
    hold on
    set(gca,'Ydir','reverse')
    if any(Zfr)
        xl = xlim;
        plot(xl,[Zfr Zfr],'r--')
    end
    legend('\sigma_{r}=p','\sigma_{z}','\sigma_{\theta}');
    xlabel('MPa')
    ylabel('Depth (km)')
    
    subplot(1,2,2);
    plot(Smax/1e6,z,Sfail/1e6,z);
    set(gca,'Ydir','reverse')
    Ds = Sfail-Smax; %<0 when failure
    diff = abs(Ds)-Ds; %will be >0 when failure, 0 elsewhere
    
    if any(diff)
        hold on;
        i1 = find(diff, 1, 'first');
        i2 = find(diff, 1, 'last');
        xl = xlim;
        plot(xl,[z(i1) z(i1)],'b--')
        plot(xl,[z(i2) z(i2)],'b--')
    end
    
    ds = Smax-Szz + (A.depth-z).*1000*A.g; %pore pressure gradient;
    match = (ds<1e-5); %finds when Smax is Szz
%     if any(match)
%         hold on;
%         plot(zeros(size(z(match))),z(match),'g*')
%         
%     end
    
    
    hold on;
    if any(Zfr)
        xl = xlim;
        plot(xl,[Zfr Zfr],'r--')
    end 
    legend('\sigma_1','Failure Criterion')
    xlabel('MPa')
    ylabel('Depth (km)')
    set(gcf,'Units','inches',...
 'Position',[0 0 9 11])
subplot(1,2,1);
xlim([-10, 200])
subplot(1,2,2);
xlim([-10, 250])
    %saveas(gcf,['op' num2str(op/1e6) '.png'])
end
