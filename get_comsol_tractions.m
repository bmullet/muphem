%% Get comsol tractions
% Last edited: 4 March 2021

% Get model
A = Amodels.initA_paper_MSH_E;
label = 'E';
chamber_fac_vec = [0.7:.02:1];

for i = 1:length(chamber_fac_vec)
   disp(["Working on ", num2str(i)])
   chamber_fac = chamber_fac_vec(i);
   A.chamber_fac = chamber_fac;
   A = Amodels.initA_paper_MSH_E(A);
   
   out = muphem("multiflow2",0,A);
   
   % Save .mat file
   save(sprintf(['./COMSOL_input/p_MSH_' label '_ch_%2d.mat'], chamber_fac*100), 'out')
   
   % Get data for comsol
 
   % srz and p
   mat = [out{2} out{17}];
   writematrix(mat, sprintf(['./COMSOL_input/p_MSH_' label '_ch_%2d.csv'], chamber_fac*100))
   mat = [out{2} out{18}];
   writematrix(mat, sprintf(['./COMSOL_input/srz_MSH_' label '_ch_%2d.csv'], chamber_fac*100))
   
   % cohesion
   mat = [out{2} A.mc.C(out{2})];
   writematrix(mat, sprintf(['./COMSOL_input/c_MSH_' label '_ch_%2d.csv'], chamber_fac*100))
   
end

%% Get fragmentation depths
label = 'E';
chamber_fac_vec = [0.7:.02:1];
frags = zeros(1,length(chamber_fac_vec));

for i = 1:length(chamber_fac_vec)
   chamber_fac = chamber_fac_vec(i);
   
   % Save .mat file
   load(sprintf(['./COMSOL_input/p_MSH_' label '_ch_%2d.mat'], chamber_fac*100))
   
   Srr = out{17};
   Srz = out{18};
   Srz = Srz;
   zvec = out{2};
   
   [m, j] = max(Srz);
   
   frags(i) = zvec(j);
   
end

%% Get comsol tractions (vary radius)
% Last edited: 25 March 2021

% Get model
A = Amodels.initA_paper_MSH_D;
label = 'D';
radius_vec = [10:10:100];

for i = 1:length(radius_vec)
   
   A.r = radius_vec(i);
   A = Amodels.initA_paper_MSH_D(A);
   
   out = muphem("multiflow2",0,A);
   
   % Save .mat file
   save(sprintf(['./COMSOL_input/p_MSH_' label '_r_%2d.mat'], radius_vec(i)), 'out')
   
   % Get data for comsol
 
   % srz and p
   mat = [out{2} out{17}];
   writematrix(mat, sprintf(['./COMSOL_input/p_MSH_' label '_r_%2d.csv'], radius_vec(i)))
   mat = [out{2} out{18}];
   writematrix(mat, sprintf(['./COMSOL_input/srz_MSH_' label '_r_%2d.csv'], radius_vec(i)))
   
   % cohesion
   mat = [out{2} A.mc.C(out{2})];
   writematrix(mat, sprintf(['./COMSOL_input/c_MSH_' label '_r_%2d.csv'], radius_vec(i)))
   
end

%% Get fragmentation depths
label = 'D';
radius_vec = [10:10:100];
frags = zeros(1,length(radius_vec));

for i = 1:length(radius_vec)
   r = radius_vec(i);
   
   % Save .mat file
   load(sprintf(['./COMSOL_input/p_MSH_' label '_r_%2d.mat'], radius_vec(i)))
   
   Srr = out{17};
   Srz = out{18};
   Srz = Srz;
   zvec = out{2};
   
   [m, j] = max(Srz);
   
   frags(i) = zvec(j);
   
end