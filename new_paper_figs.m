%%%%%%%%%%%%%%%%%%%%%
%%% Paper Figures %%%
%%%%%%%%%%%%%%%%%%%%%


%% Set default figure paarmeters
% Source: https://dgleich.wordpress.com/2013/06/04/creating-high-quality-graphics-in-matlab-for-papers-and-presentations/

% Defaults for this blog post
width = 4;     % Width in inches
height = 4;    % Height in inches
alw = 1.5;    % AxesLineWidth
fsz = 15;      % Fontsize
lw = 4;      % LineWidth
msz = 8;       % MarkerSize

% The new defaults will not take effect if there are any open figures. To
% use them, we close all figures, and then repeat the first example.
close all;

% The properties we've been using in the figures
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz
set(0,'defaultAxesFontSize',fsz)

% Set the default Size for display
defpos = get(0,'defaultFigurePosition');
set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*100]);

% Set the defaults for saving/printing to a file
set(0,'defaultFigureInvertHardcopy','on'); % This is the default anyway
set(0,'defaultFigurePaperUnits','inches'); % This is the default anyway
defsize = get(gcf, 'PaperSize');
left = (defsize(1)- width)/2;
bottom = (defsize(2)- height)/2;
defsize = [left, bottom, width, height];
set(0, 'defaultFigurePaperPosition', defsize);
set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultTextFontname', 'Arial')
%% Fragmentation depth tractions
labels = {'A', 'B'};
chamber_fac_vec = [0.7:.02:1];
frags = zeros(length(labels),length(chamber_fac_vec));
maxsrz = zeros(length(labels),length(chamber_fac_vec));
pfrag = zeros(length(labels),length(chamber_fac_vec));
A = Amodels.initA_paper_MSH_A;

for i = 1:length(labels)
    label = labels{i}
    for j = 1:length(chamber_fac_vec)
        chamber_fac = chamber_fac_vec(j);
        
        % Save .mat file
        load(sprintf(['./COMSOL_input/p_MSH_' label '_ch_%2d.mat'], chamber_fac*100))
        
        Srr = out{17};
        Srz = out{18};
        zvec = out{2};
        Szz = A.k.rho*9.806*abs(zvec);
        
        
        [m, k] = max(Srz);
        
        frags(i,j) = zvec(k);
        maxsrz(i,j) = m/Szz(k);
        pfrag(i,j) = Srr(k)./Szz(k);
        
    end
end
figure
subplot(1,3,1)
plot(chamber_fac_vec, frags(1,:)); hold on;
plot(chamber_fac_vec, frags(2,:));
legend("\phi_f = 0.80", "\phi_f = 0.70", "Location", "NorthWest")
xlabel("p_{ch}/S_{z}(z_{ch})")
ylabel("Fragmentation Depth (m)")
text(0.025,0.95,'(a)','Units','normalized','FontSize',20)
grid on

subplot(1,3,2)
plot(chamber_fac_vec, pfrag(1,:)); hold on;
plot(chamber_fac_vec, pfrag(2,:));
legend("\phi_f = 0.80", "\phi_f = 0.70", "Location", "NorthWest")
xlabel("p_{ch}/S_{z}(z_{ch})")
ylabel("p/S_z at fragmentation")
text(0.025,0.95,'(b)','Units','normalized','FontSize',20)
grid on

subplot(1,3,3)
plot(chamber_fac_vec, maxsrz(1,:)); hold on;
plot(chamber_fac_vec, maxsrz(2,:));
legend("\phi_f = 0.80", "\phi_f = 0.70", "Location", "NorthWest")
xlabel("p_{ch}/S_{z}(z_{ch})")
ylabel("\tau_{rz}/S_z at fragmentation")
text(0.025,0.95,'(c)','Units','normalized','FontSize',20)
grid on

%set(gcf, 'PaperPosition', [0 0 13 4]); %Position plot at left hand corner with width 5 and height 5.
%set(gcf, 'PaperSize', [13 4]); %Set the paper to have width 5 and height 5.
%saveas(gcf, './Figures/conduit_tractions', 'pdf') %Save figure

%% Analytical comparison at the fragmentation depth
frag_depths = 1000*[-1.4243 ,  -1.3589,   -1.2940 ,  -1.2296 ,  -1.1662 ,  -1.1035 ,  -1.0415  , -0.9803 ,  -0.9202 ,  -0.8608  , -0.8024 ,  -0.7449 ,  -0.6885 ,  -0.6333 ,  -0.5791 ,  -0.5262];
COHESION = 10e6; % MPa
A = Amodels.initA_paper_MSH_A;
label = 'A';
chamber_fac_vec = [0.7:.02:1];
plot_chamber_fac_idx =[1,6,11];

k = 0.74; % needs to be set lower too

comsol_data = importdata("./COMSOL_output/A_mc_full_wall_k_74.txt");

% prune data
ch_com = cell(length(chamber_fac_vec),1);

idx = find(diff(comsol_data(:,1)) < 0);
idx = [1; idx];
idx = [idx; length(comsol_data)];

for i = 1:length(ch_com)
    ch_com{i} = comsol_data(idx(i)+1:idx(i+1),:);
end

figure
subplot(222)
title("\phi_{frag} = 0.80"); hold on;
handles = [];
for i=plot_chamber_fac_idx
    ch = ch_com{i};
    mci = ch(1:end,2) - (15e6-COHESION)*cosd(35); %need to subtract this because in comsol we used 15MPa cohesion
    h = plot(-mci, ch(1:end,1)-5000); hold on;
    handles = [handles h];
    
end

subplot(224)
title("\phi_{frag} = 0.80"); hold on;
handles = [];
for i=plot_chamber_fac_idx
    ch = ch_com{i};
    mci = ch(1:end,2) - (15e6-COHESION)*cosd(35); %need to subtract this because in comsol we used 15MPa cohesion
    h = plot(-mci, ch(1:end,1)-5000); hold on;
    handles = [handles h];
    
end


min_an = zeros(size(chamber_fac_vec));

for j = plot_chamber_fac_idx
    chamber_fac = chamber_fac_vec(j);
    
    % Load .mat file
    load(sprintf(['./COMSOL_input/p_MSH_' label '_ch_%2d.mat'], chamber_fac*100))
    
    
    Srr = out{17};
    Srz = out{18};
    Srz = Srz;
    zvec = out{2};
    
    A = out{1};
    Szz = A.k.rho*9.806*abs(zvec);
    phi = 35/180*pi;
    pp = abs(zvec)*9.806*1000;
    %cohesion = A.mc.C(zvec) + 2e6; % COMPARISON TO COMSOL ADDS 2MPa
    cohesion = COHESION;
    
    SR = 1/2*(Srr + Szz) + ((1/2*(Srr - Szz)).^2 + Srz.^2).^(1/2);
    SZ = 1/2*(Srr + Szz) - ((1/2*(Srr - Szz)).^2 + Srz.^2).^(1/2);
    
    Stheta = 2*k.*Szz - Srr;
    
    S1 = max(max(SR,SZ), Stheta);
    S3 = min(min(SR,SZ), Stheta);
    mc = (S3-S1)/2 + (S1+S3 -2*pp)/2*sin(phi) + cohesion.*cos(phi);
    
    subplot(222)
    plot(-mc, zvec, '-k', "LineWidth", 1);
    subplot(224)
    plot(-mc, zvec, '-k', "LineWidth", 1);
    
    
    % get min around fragmentation depth
    fd = frag_depths(j);
    fdidx = (zvec > (fd - 50)) & (zvec < (fd + 50));
    min_an(j) = min(mc(fdidx));
    
    %    min_an(j) = min(mc);
    
    
end




subplot(222)
ylabel("z(m)")
plot([0, 0],ylim, '--r',  "LineWidth", 2)
xlabel("f_y(\sigma') (Pa)")

ylim([-5000, 0])
grid on;

subplot(224)
ylabel("z(m)")
plot([0, 0],ylim, '--r',  "LineWidth", 2)
xlabel("f_y(\sigma') (Pa)")

ylim(1e3*[-1.4750, -0.7510])
xlim(1e6*[-5.7216, 5.4715])
grid on;



% annotate a grid
subplot(224)
yl = ylim;
xl = xlim;


x1=xl(1);
x2=xl(2);
y1=yl(1);
y2=yl(2);
x = [x1, x2, x2, x1, x1];
y = [y1, y1, y2, y2, y1];
subplot(222)
plot(x, y, 'g--', 'LineWidth', 1);
hold on;
text(0.025,0.95,'(b)','Units','normalized','FontSize',12)


subplot(224)
plot(x, y, 'g--', 'LineWidth', 1);
hold on;
legend(handles, "0.7","0.8", "0.9");

[hleg,icons,plots] = legend("0.7", "0.8", "0.9","Location","SouthWest");
title(hleg,'p_{ch}/S_{z}(z_{ch})')
hleg.Title.Visible = 'on';
% the addition in height needed for the title:
title_hight = hleg.Position(4)/numel(plots);
hleg.Position([2 4]) = [hleg.Position(2)-title_hight hleg.Position(4)+title_hight];
hleg.Position([3]) = hleg.Position([3])*1.5;
% calculate new position for the elements in the legeng:
new_pos = fliplr(0.5/(numel(plots)+1):1/(numel(plots)+1):1);
hleg.Title.NodeChildren.Position = [0.5 new_pos(1) 0];
% set the text to the right position:
leg_txt = findobj(icons,'Type','Text');
txt_pos = cell2mat({leg_txt.Position}.');
txt_pos(:,2) = new_pos(2:end);
set(leg_txt,{'Position'},mat2cell(txt_pos,ones(numel(plots),1),3));
% set the icons to the right position:
leg_att = findobj(icons,'Type','Line');
% set(leg_att,{'YData'},mat2cell(repmat(repelem(new_pos(2:end).',...
%     numel(plots)),1,2),ones(numel(plots)*2,1),2))
set(leg_att,{'YData'},mat2cell(repmat(repelem(new_pos(2:end).',...
    2),1,2),ones(numel(plots)*2,1),2))



text(0.025,0.95,'(d)','Units','normalized','FontSize',12)
% B

frag_depths = [-1661.3, -1596.5, -1534.1, -1470.9, -1409.3, -1348.1, -1288,  -1228.3, -1169.3, -1111,  -1053.5,  -997.3,  -941.6,  -886.5,  -832.5,  -779.4]; % model B

COHESION = 10e6; % MPa


A = Amodels.initA_paper_MSH_B;
label = 'B';
chamber_fac_vec = [0.7:.02:1];
plot_chamber_fac_idx =[1,6,11];
k = 0.74; % needs to be set lower too

comsol_data = importdata("./COMSOL_output/B_mc_full_wall_k_74.txt");


% prune data
ch_com = cell(length(chamber_fac_vec),1);

idx = find(diff(comsol_data(:,1)) < 0);
idx = [1; idx];
idx = [idx; length(comsol_data)];

for i = 1:length(ch_com)
    ch_com{i} = comsol_data(idx(i)+1:idx(i+1),:);
end

subplot(221)
title("\phi_{frag} = 0.70"); hold on;
handles = [];
for i=plot_chamber_fac_idx
    ch = ch_com{i};
    mci = ch(1:end,2) - (15e6-COHESION)*cosd(35); %need to subtract this because in comsol we used 15MPa cohesion
    h = plot(-mci, ch(1:end,1)-5000); hold on;
    handles = [handles h];
end
grid on

subplot(223)
title("\phi_{frag} = 0.70"); hold on;
handles = [];
for i=plot_chamber_fac_idx
    ch = ch_com{i};
    mci = ch(1:end,2) - (15e6-COHESION)*cosd(35); %need to subtract this because in comsol we used 15MPa cohesion
    h = plot(-mci, ch(1:end,1)-5000); hold on;
    handles = [handles h];
    
end

ylim([-1710, -1000])


min_an = zeros(size(chamber_fac_vec));

for j = plot_chamber_fac_idx
    chamber_fac = chamber_fac_vec(j);
    
    % Load .mat file
    load(sprintf(['./COMSOL_input/p_MSH_' label '_ch_%2d.mat'], chamber_fac*100))
    
    Srr = out{17};
    Srz = out{18};
    Srz = Srz;
    zvec = out{2};
    
    A = out{1};
    Szz = A.k.rho*9.806*abs(zvec);
    phi = 35/180*pi;
    pp = abs(zvec)*9.806*1000;
    cohesion = COHESION;
    
    SR = 1/2*(Srr + Szz) + ((1/2*(Srr - Szz)).^2 + Srz.^2).^(1/2);
    SZ = 1/2*(Srr + Szz) - ((1/2*(Srr - Szz)).^2 + Srz.^2).^(1/2);
    
    Stheta = 2*k.*Szz - Srr;
    
    S1 = max(max(SR,SZ), Stheta);
    S3 = min(min(SR,SZ), Stheta);
    mc = (S3-S1)/2 + (S1+S3 -2*pp)/2*sin(phi) + cohesion.*cos(phi);
    
    subplot(221)
    plot(-mc, zvec, '-k', "LineWidth", 1); hold on;
    subplot(223)
    plot(-mc, zvec, '-k', "LineWidth", 1); hold on;
    
    % get min around fragmentation depth
    fd = frag_depths(j);
    fdidx = (zvec > (fd - 50)) & (zvec < (fd + 50));
    min_an(j) = min(mc(fdidx));
    
    %    min_an(j) = min(mc);
    
    
end

plot([0, 0],ylim, '--r', "LineWidth", 2)

%legend("0.70","0.72","0.74","0.76","0.78","0.80","0.82","0.84","0.86","0.88","0.90","0.92","0.94","0.96","0.98","1.00")
ylabel("z(m)")
xlabel("f_y(\sigma') (Pa)")
grid on;

subplot(221)
ylabel("z(m)")

plot([0, 0],ylim, '--r', "LineWidth", 2)

xlabel("f_y(\sigma') (Pa)")
ylim([-5000, 0])
text(0.025,0.95,'(a)','Units','normalized','FontSize',12)
grid on;

% annotate a grid
subplot(223)
yl = ylim;
xl = xlim;


x1=xl(1);
x2=xl(2);
y1=yl(1);
y2=yl(2);
x = [x1, x2, x2, x1, x1];
y = [y1, y1, y2, y2, y1];
subplot(221)
plot(x, y, 'b--', 'LineWidth', 1);
hold on;


subplot(223)
plot(x, y, 'b--', 'LineWidth', 1);
hold on;

text(0.025,0.95,'(c)','Units','normalized','FontSize',12)




%% Analytical comparison with COMSOL for model A
frag_depths = 1000*[-1.4243 ,  -1.3589,   -1.2940 ,  -1.2296 ,  -1.1662 ,  -1.1035 ,  -1.0415  , -0.9803 ,  -0.9202 ,  -0.8608  , -0.8024 ,  -0.7449 ,  -0.6885 ,  -0.6333 ,  -0.5791 ,  -0.5262];
COHESION = 10e6; % MPa
A = Amodels.initA_paper_MSH_A;
label = 'A';
chamber_fac_vec = [0.7:.02:1];
k = 0.74; % needs to be set lower too

comsol_data = importdata("./COMSOL_output/A_mc_full_wall_k_74.txt");



% prune data
ch_com = cell(length(chamber_fac_vec),1);

idx = find(diff(comsol_data(:,1)) < 0);
idx = [1; idx];
idx = [idx; length(comsol_data)];

for i = 1:length(ch_com)
    ch_com{i} = comsol_data(idx(i)+1:idx(i+1),:);
end

figure
for i=1:length(ch_com)
    ch = ch_com{i};
    mci = ch(1:end,2) - (15e6-COHESION)*cosd(35); %need to subtract this because in comsol we used 15MPa cohesion
    plot(ch(1:end,1)-5000, mci); hold on;
    
end

min_an_with_shear = zeros(size(chamber_fac_vec));
min_an_no_shear = zeros(size(chamber_fac_vec));

for j = 1:length(chamber_fac_vec)
    chamber_fac = chamber_fac_vec(j);
    
    % Load .mat file
    load(sprintf(['./COMSOL_input/p_MSH_' label '_ch_%2d.mat'], chamber_fac*100))
    
    
    Srr = out{17};
    Srz = out{18};
    Srz = Srz;
    zvec = out{2};
    
    A = out{1};
    Szz = A.k.rho*9.806*abs(zvec);
    phi = 35/180*pi;
    pp = abs(zvec)*9.806*1000;
    %cohesion = A.mc.C(zvec) + 2e6; % COMPARISON TO COMSOL ADDS 2MPa
    cohesion = COHESION;
    
    SR = 1/2*(Srr + Szz) + ((1/2*(Srr - Szz)).^2 + Srz.^2).^(1/2);
    SZ = 1/2*(Srr + Szz) - ((1/2*(Srr - Szz)).^2 + Srz.^2).^(1/2);
    
    Stheta = 2*k.*Szz - Srr;
    
    S1 = max(max(SR,SZ), Stheta);
    S3 = min(min(SR,SZ), Stheta);
    mc = (S3-S1)/2 + (S1+S3 -2*pp)/2*sin(phi) + cohesion.*cos(phi);
    
    plot(zvec, mc, '-k', "LineWidth", 1);
    
    % get min around fragmentation depth
    fd = frag_depths(j);
    fdidx = (zvec > (fd - 50)) & (zvec < (fd + 50));
    min_an_with_shear(j) = min(mc(fdidx));
    
    
    S1 = max(max(Srr,Szz), Stheta);
    S3 = min(min(Srr,Szz), Stheta);
    mc = (S3-S1)/2 + (S1+S3 -2*pp)/2*sin(phi) + cohesion.*cos(phi);
    
    
    % get min around fragmentation depth
    fd = frag_depths(j);
    fdidx = (zvec > (fd - 50)) & (zvec < (fd + 50));
    min_an_no_shear(j) = min(mc(fdidx));
    
    %    min_an(j) = min(mc);
    
    
end

plot(xlim,[0, 0], '--r')

xlabel("z(m)")
ylabel("mc (Pa)")

%% Minimums
min_com = zeros(size(ch_com));
for i = 1:length(ch_com)
    ch = ch_com{i};
    fd = frag_depths(i);
    fdidx = (ch(:,1) > (fd - 50 + 5000)) & (ch(:,1) < (fd + 50 + 5000));
    min_com(i) = min(ch(fdidx,2)) - (15e6 - COHESION)*cos(phi);
    %min_com(i) = min(ch(:,2));
end

figure(5); subplot(122);
h1 = plot(chamber_fac_vec, -min_com, "*k"); hold on;
h2 = plot(chamber_fac_vec, -min_an_with_shear);
h3 = plot(chamber_fac_vec, -min_an_no_shear);

xlabel("p_{ch}/S_z(z_{ch})")
ylabel("Max(f_y(\sigma') (Pa)")
plot(xlim, [0,0], '--r')

legend([h3, h2, h1], "Analytical - no shear", "Analytical - with shear", "FEM");
title("\phi_{frag} = 0.80")
text(0.9,0.95,'(b)','Units','normalized','FontSize',20)
grid on;



%% Compare across k
label = 'A';
comsol_data = importdata("./COMSOL_output/A_mc_full_wall_all_k.txt");
frag_depths = 1000*[-1.4243 ,  -1.3589,   -1.2940 ,  -1.2296 ,  -1.1662 ,  -1.1035 ,  -1.0415  , -0.9803 ,  -0.9202 ,  -0.8608  , -0.8024 ,  -0.7449 ,  -0.6885 ,  -0.6333 ,  -0.5791 ,  -0.5262];
COHESION = 10e6; % MPa
chamber_fac_vec = [0.7, 0.8, 0.9];

idx = [1,6,11];
frag_depths = frag_depths(idx);

colors = lines(3);

ch_com = cell(length(idx),1);

for i=1:length(chamber_fac_vec)
    ch_com{i} = comsol_data(find(comsol_data(:,2) == chamber_fac_vec(i)),:);
end

figure(6);
subplot(122)
for i=1:length(ch_com)
    ch = ch_com{i};
    plot(ch(1:end,1), -(ch(1:end,3) - (15e6-COHESION)*cosd(35))); hold on;
    chamber_fac_vec(i) = ch(1,2);
end

kk = 0.50:.01:1;


for j = 1:length(chamber_fac_vec)
    mcvec = zeros(size(kk));
    chamber_fac = chamber_fac_vec(j);
    
    % Load .mat file
    load(sprintf(['./COMSOL_input/p_MSH_' label '_ch_%2d.mat'], chamber_fac*100))
    
    Srr = out{17};
    Srz = out{18};
    zvec = out{2};
    
    A = out{1};
    Szz = A.k.rho*9.806*abs(zvec);
    phi = 35/180*pi;
    pp = abs(zvec)*9.806*1000;
    cohesion = COHESION; % COMPARISON TO COMSOL ADDS 2MPa
    
    SR = 1/2*(Srr + Szz) + ((1/2*(Srr - Szz)).^2 + Srz.^2).^(1/2);
    SZ = 1/2*(Srr + Szz) - ((1/2*(Srr - Szz)).^2 + Srz.^2).^(1/2);
    for l = 1:length(kk)
        k = kk(l);
        Stheta = 2*k.*Szz - Srr;
        
        S1 = max(max(SR,SZ), Stheta);
        S3 = min(min(SR,SZ), Stheta);
        
        mc = (S3-S1)/2 + (S1+S3 -2*pp)/2*sin(phi) + cohesion*cos(phi);
        fd = frag_depths(j);
        fdidx = (zvec > (fd - 50)) & (zvec < (fd + 50));
        mcvec(l) = min(mc(fdidx));
        %mcvec(l) = min(mc);
        
    end
    
    plot(kk, -mcvec, "Color", colors(j,:), "LineStyle", "--")
    
end

xlabel('k')
ylabel("Max(f_y(\sigma') (Pa)")

grid on;
ylim([-5e6,10e6])
xlim([0.5, 1])
title("\phi_{frag} = 0.8")

[hleg,icons,plots] = legend("0.7", "0.8", "0.9","Location","NorthWest");
title(hleg,'p_{ch}/S_{z}(z_{ch})')
hleg.Title.Visible = 'on';
% the addition in height needed for the title:
title_hight = hleg.Position(4)/numel(plots);
hleg.Position([2 4]) = [hleg.Position(2)-title_hight hleg.Position(4)+title_hight];
hleg.Position([3]) = hleg.Position([3])*1.5;
% calculate new position for the elements in the legeng:
new_pos = fliplr(0.5/(numel(plots)+1):1/(numel(plots)+1):1);
hleg.Title.NodeChildren.Position = [0.5 new_pos(1) 0];
% set the text to the right position:
leg_txt = findobj(icons,'Type','Text');
txt_pos = cell2mat({leg_txt.Position}.');
txt_pos(:,2) = new_pos(2:end);
set(leg_txt,{'Position'},mat2cell(txt_pos,ones(numel(plots),1),3));
% set the icons to the right position:
leg_att = findobj(icons,'Type','Line');
% set(leg_att,{'YData'},mat2cell(repmat(repelem(new_pos(2:end).',...
%     numel(plots)),1,2),ones(numel(plots)*2,1),2))
set(leg_att,{'YData'},mat2cell(repmat(repelem(new_pos(2:end).',...
    2),1,2),ones(numel(plots)*2,1),2))

text(0.025,0.95,'(b)','Units','normalized','FontSize',20)


%% Analytical comparison with COMSOL for model B

frag_depths = [-1661.3, -1596.5, -1534.1, -1470.9, -1409.3, -1348.1, -1288,  -1228.3, -1169.3, -1111,  -1053.5,  -997.3,  -941.6,  -886.5,  -832.5,  -779.4]; % model B

A = Amodels.initA_paper_MSH_B;
label = 'B';
chamber_fac_vec = [0.7:.02:1];
k = 0.70; % needs to be set lower too

comsol_data = importdata("./COMSOL_output/B_mc_full_wall_k_70_3.txt");


% prune data
ch_com = cell(length(chamber_fac_vec),1);

idx = find(diff(comsol_data(:,1)) < 0);
idx = [1; idx];
idx = [idx; length(comsol_data)];

for i = 1:length(ch_com)
    ch_com{i} = comsol_data(idx(i)+1:idx(i+1),:);
end

figure
for i=1:length(ch_com)
    ch = ch_com{i};
    plot(ch(1:end,1)-5000, ch(1:end,2) - (15e6 - COHESION)*cos(phi)); hold on;
    
end

min_an_with_shear = zeros(size(chamber_fac_vec));
min_an_no_shear = zeros(size(chamber_fac_vec));

for j = 1:length(chamber_fac_vec)
    chamber_fac = chamber_fac_vec(j);
    
    % Load .mat file
    load(sprintf(['./COMSOL_input/p_MSH_' label '_ch_%2d.mat'], chamber_fac*100))
    
    Srr = out{17};
    Srz = out{18};
    Srz = Srz;
    zvec = out{2};
    
    A = out{1};
    Szz = A.k.rho*9.806*abs(zvec);
    phi = 35/180*pi;
    pp = abs(zvec)*9.806*1000;
    cohesion = COHESION;
    
    SR = 1/2*(Srr + Szz) + ((1/2*(Srr - Szz)).^2 + Srz.^2).^(1/2);
    SZ = 1/2*(Srr + Szz) - ((1/2*(Srr - Szz)).^2 + Srz.^2).^(1/2);
    
    Stheta = 2*k.*Szz - Srr;
    
    S1 = max(max(SR,SZ), Stheta);
    S3 = min(min(SR,SZ), Stheta);
    mc = (S3-S1)/2 + (S1+S3 -2*pp)/2*sin(phi) + cohesion.*cos(phi);
    
    plot(zvec, mc, '-k', "LineWidth", 1); hold on;
    
    % get min around fragmentation depth
    fd = frag_depths(j);
    fdidx = (zvec > (fd - 50)) & (zvec < (fd + 50));
    min_an_with_shear(j) = min(mc(fdidx));
    
    
    S1 = max(max(Srr,Szz), Stheta);
    S3 = min(min(Srr,Szz), Stheta);
    mc = (S3-S1)/2 + (S1+S3 -2*pp)/2*sin(phi) + cohesion.*cos(phi);
    
    % get min around fragmentation depth
    fd = frag_depths(j);
    fdidx = (zvec > (fd - 50)) & (zvec < (fd + 50));
    min_an_no_shear(j) = min(mc(fdidx));
    
    
    %    min_an(j) = min(mc);
    
    
end

plot(xlim,[0, 0], '--r')

legend("0.70","0.72","0.74","0.76","0.78","0.80","0.82","0.84","0.86","0.88","0.90","0.92","0.94","0.96","0.98","1.00")
xlabel("z(m)")
ylabel("mc (Pa)")

%% Minimums
frag_depths = [-1661.3, -1596.5, -1534.1, -1470.9, -1409.3, -1348.1, -1288,  -1228.3, -1169.3, -1111,  -1053.5,  -997.3,  -941.6,  -886.5,  -832.5,  -779.4]; % model B

min_com = zeros(size(ch_com));

for i = 1:length(ch_com)
    ch = ch_com{i};
    fd = frag_depths(i);
    fdidx = (ch(:,1) > (fd - 50 + 5000)) & (ch(:,1) < (fd + 50 + 5000));
    min_com(i) = min(ch(fdidx,2)) - (15e6-COHESION)*cos(phi);
    %min_com(i) = min(ch(:,2));
end

figure(5); subplot(121);
h1 = plot(chamber_fac_vec, -min_com, "*k"); hold on;
h2 = plot(chamber_fac_vec, -min_an_with_shear);
h3 = plot(chamber_fac_vec, -min_an_no_shear);

xlabel("p_{ch}/S_z(z_{ch})")
ylabel("Max(f_y(\sigma') (Pa)")
plot(xlim, [0,0], '--r')

legend([h3, h2, h1], "Analytical - no shear", "Analytical - with shear", "FEM");
title("\phi_{frag} = 0.70")
grid on;

text(0.9,0.95,'(a)','Units','normalized','FontSize',20)

% set(gcf, 'PaperPosition', [0 0 9 4]); %Position plot at left hand corner with width 5 and height 5.
% set(gcf, 'PaperSize', [9 4]); %Set the paper to have width 5 and height 5.
% saveas(gcf, './Figures/min_mc_over_pch', 'pdf') %Save figure

%% Compare across k
label = 'B';
comsol_data = importdata("./COMSOL_output/B_mc_frag_depth_all_k.txt");
frag_depths = [-1661.3, -1596.5, -1534.1, -1470.9, -1409.3, -1348.1, -1288,  -1228.3, -1169.3, -1111,  -1053.5,  -997.3,  -941.6,  -886.5,  -832.5,  -779.4]; % model B

COHESION = 10e6; % MPa
chamber_fac_vec = [0.7, 0.8, 0.9];

idx = [1,6,11];
frag_depths = frag_depths(idx);

colors = lines(3);

ch_com = cell(length(idx),1);

for i=1:length(chamber_fac_vec)
    ch_com{i} = comsol_data(find(comsol_data(:,2) == chamber_fac_vec(i)),:);
end


figure(6)
subplot(121)
title('\phi_{frag} = 0.7'); hold on;
handles = [];
for i=1:length(ch_com)
    ch = ch_com{i};
    h = plot(ch(1:end,1), -(ch(1:end,3) - (15e6-COHESION)*cosd(35))); hold on;
    handles = [handles h];
    chamber_fac_vec(i) = ch(1,2);
end

kk = 0.50:.01:1;

for j = 1:length(chamber_fac_vec)
    mcvec = zeros(size(kk));
    chamber_fac = chamber_fac_vec(j);
    
    % Load .mat file
    load(sprintf(['./COMSOL_input/p_MSH_' label '_ch_%2d.mat'], chamber_fac*100))
    
    Srr = out{17};
    Srz = out{18};
    zvec = out{2};
    
    A = out{1};
    Szz = A.k.rho*9.806*abs(zvec);
    phi = 35/180*pi;
    pp = abs(zvec)*9.806*1000;
    cohesion = COHESION;
    
    SR = 1/2*(Srr + Szz) + ((1/2*(Srr - Szz)).^2 + Srz.^2).^(1/2);
    SZ = 1/2*(Srr + Szz) - ((1/2*(Srr - Szz)).^2 + Srz.^2).^(1/2);
    for l = 1:length(kk)
        k = kk(l);
        Stheta = 2*k.*Szz - Srr;
        
        S1 = max(max(SR,SZ), Stheta);
        S3 = min(min(SR,SZ), Stheta);
        
        mc = (S3-S1)/2 + (S1+S3 -2*pp)/2*sin(phi) + cohesion*cos(phi);
        fd = frag_depths(j);
        fdidx = (zvec > (fd - 50)) & (zvec < (fd + 50));
        mcvec(l) = min(mc(fdidx));
        %mcvec(l) = min(mc);
        
    end
    
    plot(kk, -mcvec, "Color", colors(j,:), "LineStyle", "--")
    
end

xlabel('k')
legend(handles, "0.7", "0.8", "0.9", "Location", "SouthWest");
ylim([-5e6,10e6])
xlim([0.5,1])
ylabel("Max(f_y(\sigma') (Pa)")
grid on

text(0.025,0.95,'(a)','Units','normalized','FontSize',20)


[hleg,icons,plots] = legend("0.7", "0.8", "0.9","Location","NorthWest");
title(hleg,'p_{ch}/S_{z}(z_{ch})')
hleg.Title.Visible = 'on';
% the addition in height needed for the title:
title_hight = hleg.Position(4)/numel(plots);
hleg.Position([2 4]) = [hleg.Position(2)-title_hight hleg.Position(4)+title_hight];
hleg.Position([3]) = hleg.Position([3])*1.5;
% calculate new position for the elements in the legeng:
new_pos = fliplr(0.5/(numel(plots)+1):1/(numel(plots)+1):1);
hleg.Title.NodeChildren.Position = [0.5 new_pos(1) 0];
% set the text to the right position:
leg_txt = findobj(icons,'Type','Text');
txt_pos = cell2mat({leg_txt.Position}.');
txt_pos(:,2) = new_pos(2:end);
set(leg_txt,{'Position'},mat2cell(txt_pos,ones(numel(plots),1),3));
% set the icons to the right position:
leg_att = findobj(icons,'Type','Line');
% set(leg_att,{'YData'},mat2cell(repmat(repelem(new_pos(2:end).',...
%     numel(plots)),1,2),ones(numel(plots)*2,1),2))
set(leg_att,{'YData'},mat2cell(repmat(repelem(new_pos(2:end).',...
    2),1,2),ones(numel(plots)*2,1),2))


% set(gcf, 'PaperPosition', [0 0 9 4]); %Position plot at left hand corner with width 5 and height 5.
% set(gcf, 'PaperSize', [9 4]); %Set the paper to have width 5 and height 5.
% saveas(gcf, './Figures/frag_failure_all_k', 'pdf') %Save figure

%% Look at matlab solns

figure
chamber_fac_vec_const = [0.7:.02:1];
colors = lines(length(chamber_fac_vec_const));
handles = [];

for j = 1:length(chamber_fac_vec_const)
    %mcvec = zeros(size(kk));
    chamber_fac = chamber_fac_vec_const(j);
    
    % Load .mat file
    load(sprintf(['./COMSOL_input/p_MSH_' label '_ch_%2d.mat'], chamber_fac*100))
    
    Srr = out{17};
    Srz = out{18};
    zvec = out{2};
    
    plot(Srr, zvec, "Color", colors(j,:)); hold on;
    h = plot(Srz, zvec, "Color", colors(j,:)); hold on;
    handles = [handles h];
    
    %     A = out{1};
    %     Szz = A.k.rho*9.806*abs(zvec);
    %     phi = 35/180*pi;
    %     pp = abs(zvec)*9.806*1000;
    %     cohesion = A.mc.C(zvec) + 2e6; % COMPARISON TO COMSOL ADDS 2MPa
    %
    %     SR = 1/2*(Srr + Szz) + ((1/2*(Srr - Szz)).^2 + Srz.^2).^(1/2);
    %     SZ = 1/2*(Srr + Szz) - ((1/2*(Srr - Szz)).^2 + Srz.^2).^(1/2);
    %     for l = 1:length(kk)
    %         k = kk(l);
    %         Stheta = 2*k.*Szz - Srr;
    %
    %         S1 = max(max(SR,SZ), Stheta);
    %         S3 = min(min(SR,SZ), Stheta);
    %
    %         mc = (S3-S1)/2 + (S1+S3 -2*pp)/2*sin(phi) + cohesion*cos(phi);
    %
    %         mcvec(l) = min(mc);
    %
    %     end
    
    %    plot(kk, mcvec, "Color", colors(j,:), "LineStyle", "--")
    
end
legend(handles, "0.70","0.72","0.74","0.76","0.78","0.80","0.82","0.84","0.86","0.88","0.90","0.92","0.94","0.96","0.98","1.00")

%% Vary radius

% Load comsol data

A = Amodels.initA_paper_MSH_D;
label = 'D';
radius_vec = [10:10:100];
frags = zeros(1,length(radius_vec));

comsol_data = importdata("./COMSOL_output/D_mc_full_wall_vary_r.txt");

% prune data
ch_com = cell(length(radius_vec),1);

idx = find(diff(comsol_data(:,1)) < 0);
idx = [1; idx];
idx = [idx; length(comsol_data)];

for i = 1:length(ch_com)
    ch_com{i} = comsol_data(idx(i)+1:idx(i+1),:);
end

figure
colors = parula(length(ch_com));
for i=1:length(ch_com)
    ch = ch_com{i};
    plot( -(ch(1:end,2) - (15e6 - COHESION)*cos(phi)), ch(1:end,1)-5000); hold on;
    
end

k = 0.74;

radius_vec = [10:10:100];
frags = zeros(1,length(radius_vec));
rzp = zeros(1, length(radius_vec));
tauplot = zeros(1, length(radius_vec));
pplot = zeros(1, length(radius_vec));
tauzplot = zeros(1, length(radius_vec));
mc_an = zeros(1, length(radius_vec));
mc_com = zeros(1, length(radius_vec));

for i = 1:length(radius_vec)
    r = radius_vec(i);
    
    % Save .mat file
    load(sprintf(['./COMSOL_input/p_MSH_' label '_r_%2d.mat'], radius_vec(i)))
    Srr = out{17};
    Srz = out{18};
    zvec = out{2};
    Szz = A.k.rho*9.806*abs(zvec);
    
    [rzp(i), l] = max(Srz./Srr);
    
    tauplot(i) = Srz(l);
    pplot(i) = Srr(l);
    tauzplot(i) = Srz(l)./Szz(l);
    
    A = out{1};
    phi = 35/180*pi;
    pp = abs(zvec)*9.806*1000;
    cohesion = COHESION; % COMPARISON TO COMSOL ADDS 2MPa

    Stheta = 2*k.*Szz - Srr;
    
    SR = Srr;
    SZ = Szz;
    
    S1 = max(max(SR,SZ), Stheta);
    S3 = min(min(SR,SZ), Stheta);
  
    mc_an_i = -1*((S3-S1)/2 + (S1+S3 -2*pp)/2*sin(phi) + cohesion*cos(phi));
    
    ch = ch_com{i};
    mc_com_i = (-(ch(1:end,2) - (15e6 - COHESION)*cos(phi)));
    
    mc_an(i) = mc_an_i(l);
    mc_com(i) = max(mc_com_i);
    
    
end


ylabel("z(m)")
xlabel("f_y(\sigma') (Pa)")

[hleg,icons,plots] = legend("10","20","30","40","50","60","70","80","90","100","Location","SouthWest", 'FontSize',10);
hleg.FontSize = 8;
title(hleg,'R (m)')
hleg.Title.Visible = 'on';
% the addition in height needed for the title:
title_hight = hleg.Position(4)/numel(plots);
hleg.Position([2 4]) = [hleg.Position(2)-title_hight hleg.Position(4)+title_hight];
hleg.Position([3]) = hleg.Position([3])*1.5;
% calculate new position for the elements in the legeng:
new_pos = fliplr(0.5/(numel(plots)+1):1/(numel(plots)+1):1);
hleg.Title.NodeChildren.Position = [0.5 new_pos(1) 0];
% set the text to the right position:
leg_txt = findobj(icons,'Type','Text');
txt_pos = cell2mat({leg_txt.Position}.');
txt_pos(:,2) = new_pos(2:end);
set(leg_txt,{'Position'},mat2cell(txt_pos,ones(numel(plots),1),3));
% set the icons to the right position:
leg_att = findobj(icons,'Type','Line');
% set(leg_att,{'YData'},mat2cell(repmat(repelem(new_pos(2:end).',...
%     numel(plots)),1,2),ones(numel(plots)*2,1),2))
set(leg_att,{'YData'},mat2cell(repmat(repelem(new_pos(2:end).',...
    2),1,2),ones(numel(plots)*2,1),2))

plot([0,0],ylim,"--r")




figure
subplot(121)
plot(radius_vec, rzp); hold on;
plot(radius_vec, tauzplot);
xlabel("Radius (m)")
legend("\tau / p", "\tau / S_{z}")

subplot(122)
plot(radius_vec, tauplot); hold on;
plot(radius_vec, pplot);

xlabel("Radius (m)")
ylabel("Pa")
legend("\tau", "p")

figure
subplot(121)
plot(radius_vec, mc_an); hold on
plot(radius_vec, mc_com);
legend("Analytical", "FEM")
ylabel("Max(f_y(\sigma')) (Pa)")
xlabel("R (m)")
grid on
text(0.025,0.95,'(a)','Units','normalized','FontSize',20)
subplot(122); yyaxis left
plot(radius_vec, mc_com-mc_an, "*", "Color", "#0384fc")
ylabel("f_y(\sigma')_{FEM} - f_y(\sigma')_{an}")
xlabel("R (m)")
ylim([0, 3.5e6])
 yyaxis right
plot(radius_vec, tauzplot, 'Color', '#fc6703');
xlabel("R (m)")
ylabel("\tau_{rz}/S_z at fragmentation")
text(0.025,0.95,'(b)','Units','normalized','FontSize',20)

ylim([0.05, 0.25])