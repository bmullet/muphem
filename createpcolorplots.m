%% Create pcolor plots from representative data
set(0, 'DefaultAxesFontWeight', 'normal', ...
      'DefaultAxesFontSize', 20, ...
      'DefaultAxesTitleFontWeight', 'bold', ...
      'DefaultAxesTitleFontSizeMultiplier', 1, ...
      'DefaultLineLineWidth',3) ;
load('representative.mat');  
z = 1:1:3998;
pmat = zeros(length(z),41);
ummat = zeros(length(z),41);
ugmat = zeros(length(z),41);
rhogmat = zeros(length(z),41);
phimat = zeros(length(z),41);
diffmat = zeros(length(z),41);
pchvec = zeros(1,41);
Qbotvec = pchvec;

for i=30:-1:1
    out = outvec{i};
    A = out{1}; zvec = out{2}; pvec = out{3}; ugvec = out{4}; umvec = out{5};
    phivec = out{6}; rhogvec = out{7}; chidvec = out{8}; Qmvec = out{9}; Qgvec = out{10}; failure = out{11};
    Qvec = Qmvec+Qgvec;
    kill = diff(zvec)==0;
    zvec(kill)=[];
    pvec(kill)=[];
    umvec(kill) = [];
    ugvec(kill) = [];
    rhogvec(kill) = [];
    phivec(kill) = [];
    Qvec(kill)=[];

    pmat(:,i) = interp1(zvec,pvec,z);
    ummat(:,i) = interp1(zvec,umvec,z);
    ugmat(:,i) = interp1(zvec,ugvec,z);
    rhogmat(:,i) = interp1(zvec,rhogvec,z);
    phimat(:,i) = interp1(zvec,phivec,z);
    pchvec(i) = (A.Pchamber - A.depth*A.g*A.rhom0)/1e6;
    Qbotvec(i) = Qvec(1);
    
    [Srr, Szz, Stt] = kirsch(zvec,pvec,A);
    [Smax,Sfail,failure] = mcfailure(A,Srr,Szz,Stt,zvec);
    Ds = Sfail-Smax; %<0 when failure
    diffvec = abs(Ds)-Ds; %will be >0 when failure, 0 elsewhere
    if any(diffvec)
        disp('failure!')
    end
    diffmat(:,i) = interp1(zvec,diffvec,z);
    plot2dfailure(zvec,pvec,A);
    pause;
  
    disp(min(pvec))
    
end

zprint = fliplr(z)/1e3;
plot(Qbotvec,pchvec);


%% Failure region
h = pcolor(pchvec,zprint,diffmat);
set(gca,'Ydir','reverse')
set(h, 'EdgeColor', 'none');
shading interp
colormap(jet(256))
colorbar
hTitle = title('Magma Pressure (MPa)');
hXLabel = xlabel('Chamber Overpressure (MPa)  ');
hYLabel = ylabel('Depth (km)');

set([hXLabel, hYLabel]  , ...
    'FontSize'   , 20          );
set( hTitle                    , ...
    'FontSize'   , 22          , ...
    'FontWeight' , 'bold'      );

set(gcf,'Units','inches',...
 'Position',[0 0 5 7])

%% Pressure

pmat = pmat/1e6;
h = pcolor(pchvec,zprint,pmat);
set(gca,'Ydir','reverse')
set(h, 'EdgeColor', 'none');
shading interp
colormap(jet(256))
colorbar
hTitle = title('Magma Pressure (MPa)');
hXLabel = xlabel('Chamber Overpressure (MPa)  ');
hYLabel = ylabel('Depth (km)');

set([hXLabel, hYLabel]  , ...
    'FontSize'   , 20          );
set( hTitle                    , ...
    'FontSize'   , 22          , ...
    'FontWeight' , 'bold'      );

set(gcf,'Units','inches',...
 'Position',[0 0 5 7])

% %% Gas Velocity
% figure
% h = pcolor(pchvec,zprint,ugmat);
% set(gca,'Ydir','reverse')
% set(h, 'EdgeColor', 'none');
% shading interp
% colormap(jet(256))
% colorbar
% hTitle = title('Gas Velocity (m/s)');
% hXLabel = xlabel('Chamber Overpressure (MPa)  ');
% hYLabel = ylabel('Depth (km)');
% 
% set([hXLabel, hYLabel]  , ...
%     'FontSize'   , 20          );
% set( hTitle                    , ...
%     'FontSize'   , 22          , ...
%     'FontWeight' , 'bold'      );
% 
% set(gcf,'Units','inches',...
%  'Position',[0 0 5 7])
% 
% %% Melt Velocity
% figure
% h = pcolor(pchvec,zprint,ummat);
% set(gca,'Ydir','reverse')
% set(h, 'EdgeColor', 'none');
% shading interp
% colormap(jet(256))
% colorbar
% hTitle = title('Melt Velocity (m/s)');
% hXLabel = xlabel('Chamber Overpressure (MPa)  ');
% hYLabel = ylabel('Depth (km)');
% 
% set([hXLabel, hYLabel]  , ...
%     'FontSize'   , 20          );
% set( hTitle                    , ...
%     'FontSize'   , 22          , ...
%     'FontWeight' , 'bold'      );
% 
% set(gcf,'Units','inches',...
%  'Position',[0 0 5 7])
% %% Gas Density
% figure
% h = pcolor(pchvec,zprint,rhogmat);
% set(gca,'Ydir','reverse')
% set(h, 'EdgeColor', 'none');
% shading interp
% colormap(jet(256))
% colorbar
% hTitle = title('Gas Density (kg/m^3)');
% hXLabel = xlabel('Chamber Overpressure (MPa)  ');
% hYLabel = ylabel('Depth (km)');
% 
% set([hXLabel, hYLabel]  , ...
%     'FontSize'   , 20          );
% set( hTitle                    , ...
%     'FontSize'   , 22          , ...
%     'FontWeight' , 'bold'      );
% 
% set(gcf,'Units','inches',...
%  'Position',[0 0 5 7])
% 
% %% Gas Volume Fraction
% figure
% h = pcolor(pchvec,zprint,phimat);
% set(gca,'Ydir','reverse')
% set(h, 'EdgeColor', 'none');
% shading interp
% colormap(jet(256))
% colorbar
% hTitle = title('Gas Volume Fraction');
% hXLabel = xlabel('Chamber Overpressure (MPa)  ');
% hYLabel = ylabel('Depth (km)');
% 
% set([hXLabel, hYLabel]  , ...
%     'FontSize'   , 20          );
% set( hTitle                    , ...
%     'FontSize'   , 22          , ...
%     'FontWeight' , 'bold'      );
% 
% set(gcf,'Units','inches',...
%  'Position',[0 0 5 7])
