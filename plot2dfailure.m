function [Srr, Szz, Stt] = plot2dfailure (zvec,p,A)
    N = 30;
    cR = A.r;
    r = linspace(cR,cR*1.5,N);
    [Rdis,P] = meshgrid(r,p);
    
    z = A.depth - zvec; % Change into depth 
       
    Szz = (1.01e5+A.k.rho*(z)*A.g)*ones(1,N);   
    %S = Szz;            % Isotropic stress condition
    S = A.k.K*Szz;      % Poisson's Ratio stress condition
    
    Srr = S.*(1-cR^2./Rdis.^2)+P.*(cR^2./Rdis.^2); 
    Stt = (S).*(1+cR^2./Rdis.^2)-P.*(cR^2./Rdis.^2);
    ppore = (A.depth-zvec).*1000*A.g; %pore pressure gradient
    [~,PPORE] = meshgrid(r,ppore);
    
    % Find max/min
    c = nan([size(Szz),3]);
    c(:,:,1) = Szz;
    c(:,:,2) = Stt;
    c(:,:,3) = Srr;
    Smax = max(c,[],3);
    Smin = min(c,[],3);
    
    qu = 2*A.mc.C*tan(deg2rad(45)+A.mc.phi/2);
    Sfail = qu + (Smin-PPORE)*tan(deg2rad(45)+A.mc.phi/2)^2;
    Smaxef = Smax-PPORE;
    
    figure(1)
    h = pcolor(Rdis/cR,z*ones(1,N),Smaxef./Sfail);
    set(gca,'Ydir','reverse')
    set(h, 'EdgeColor', 'none');
    shading interp
    colormap('jet')
    colorbar
    caxis([0,2])
%     
%     figure(2)
%     h = pcolor(Rdis/cR,z*ones(1,N),Smaxef);
%     set(gca,'Ydir','reverse')
%     set(h, 'EdgeColor', 'none');
%     shading interp
%     colormap(jet(256))
%     colorbar
%     
%     figure(3)
%     h = pcolor(Rdis/cR,z*ones(1,N),Sfail);
%     set(gca,'Ydir','reverse')
%     set(h, 'EdgeColor', 'none');
%     shading interp
%     colormap(jet(256))
%     colorbar
%     
%     h = pcolor(Rdis/cR,z*ones(1,N),Szz-P);
%     set(gca,'Ydir','reverse')
%     set(h, 'EdgeColor', 'none');
%     shading interp
%     colormap(jet(256))
%     colorbar
%     pause
    figure(2)
    h = plot(Szz(:,1),z,'r',Stt(:,1),z,'b',Srr(:,1),z,'k');
    set(gca,'Ydir','reverse')
    legend('Sz','Stheta','Sr')
end