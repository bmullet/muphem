%load("comsol_failure.mat")

Srr = out{17};
Srz = out{18};
zvec = out{2};
A = out{1};
Szz = A.k.rho*9.806*abs(zvec);


phi = 35/180*pi;
pp = abs(zvec)*9.806*1000;
cohesion = A.mc.C(zvec); 
mu = tan(phi);

C = 2*cohesion.*((mu^2 + 1)^(1/2) + mu);


SR = 1/2*(Srr + Szz) + ((1/2*(Srr - Szz)).^2 + Srz.^2).^(1/2);
SZ = 1/2*(Srr + Szz) - ((1/2*(Srr - Szz)).^2 + Srz.^2).^(1/2);

Stheta = 2*A.lambda(zvec).*Szz - Srr;

S1 = max(max(SR,SZ), Stheta);
S3 = min(min(SR,SZ), Stheta);

% plot(zvec, Srr - Szz); hold on;
% plot(zvec, S3 - S1);
% plot(dstress_no_shear(:,1), dstress_no_shear(:,2));
% plot(dstress(:,1), dstress(:,2));
% legend('No shear', 'With shear', "COMSOL no shear", "COMSOL with shear")
% xlabel("z(m)"); ylabel("Differential stress")


%% Trim to fragmentation depth
idx = (zvec < -1200) & (zvec > -1600);
Szz = Szz(idx);
Srr = Srr(idx);
SR = SR(idx);
SZ = SZ(idx);
Srz = Srz(idx);
pp = pp(idx);
C = C(idx);
zvec = zvec(idx);

%% min mc v k
kk = 0.55:.01:1;
mcvec = zeros(size(kk));



for i = 1:length(kk)
  k = kk(i);
  Stheta = 2*k.*Szz - Srr;  
  
  S1 = max(max(SR,SZ), Stheta);
  S3 = min(min(SR,SZ), Stheta);
  
  mc = (S3-S1)/2 + (S1+S3 -2*pp)/2*sin(phi) + A.mc.C(zvec)*cos(phi);
  
  mcvec(i) = min(mc);
  
end
figure

plot(kk, mcvec)

d = importdata("./COMSOL_output/mc_failure_ch9_phi7_R50.txt");
hold on;
plot(d(:,1), d(:,2));
legend("analytical", "COMSOL")
ylabel("Mohr-Coulomb failure")
xlabel("k")

%% get val at fragdepth
[val, idx] = max(Srz);
fprintf("C/sig_z %.2f\n", C(idx)./Szz(idx))
fprintf("p/sig_z %.2f\n", Srr(idx)./Szz(idx))
fprintf("Srz/p %.2f\n", Srz(idx)./Srr(idx))


%% mc v z
load("p_ch_1_phi_60_R80.csv")
d = importdata("./COMSOL_output/mc_failure_ch1_phi6_R80.txt");

Srr = out{17};
Srz = out{18};
zvec = out{2};
A = out{1};
Szz = A.k.rho*9.806*abs(zvec);


phi = 35/180*pi;
pp = abs(zvec)*9.806*1000;
cohesion = A.mc.C(zvec); 
mu = tan(phi);

C = 2*cohesion.*((mu^2 + 1)^(1/2) + mu);


SR = 1/2*(Srr + Szz) + ((1/2*(Srr - Szz)).^2 + Srz.^2).^(1/2);
SZ = 1/2*(Srr + Szz) - ((1/2*(Srr - Szz)).^2 + Srz.^2).^(1/2);

Stheta = 2*A.lambda(zvec).*Szz - Srr;

S1 = max(max(SR,SZ), Stheta);
S3 = min(min(SR,SZ), Stheta);

mc = (S3-S1)/2 + (S1+S3 -2*pp)/2*sin(phi) + A.mc.C(zvec)*cos(phi);
  
figure

plot(zvec, mc)

hold on;
plot(d(:,1), d(:,2));
legend("analytical", "COMSOL")
ylabel("Mohr-Coulomb failure")
xlabel("k")


%%
% Get data for comsol
%load("example_eruption.mat")
% srz and p
mat = [out{2} out{17}];
writematrix(mat, 'p_ch_1_phi_60_R80.csv')
mat = [out{2} out{18}];
writematrix(mat, 'srz_ch_1_phi_60_R80.csv')

% cohesion
mat = [out{2} A.mc.C(out{2})];
writematrix(mat, 'cohesion_ch_1_phi_60_R80.csv')
