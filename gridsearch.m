% gridsearch.m
% Used to explore the solution space for the muphem model

%% CONDUCT GRID SEARCH

% Get list of json input files
json_dir = '../multiphase_test';
jsons = dir(json_dir);
jsons = {jsons.name};
idx = cellfun(@(x) contains(x,'.json'),jsons);
jsons = jsons(idx);

A = Amodels.initA_Degruyter2012;

for i=1:length(jsons)
    % Run muphem code
    input_json = ['../multiphase_test/' jsons{i}];
    params = jsondecode(fileread(input_json));
    try
        out = muphem('multiflow2',0,A,params);

        % Save output
        output_name = [input_json(1:end-5) '.mat'];
        save(output_name, 'out');
    catch
    end

end

%% SEE RESULTS
mat_dir = '../multiphase_test';
mat_files = dir(mat_dir);
mat_files = {mat_files.name};
idx = cellfun(@(x) contains(x,'.mat'),mat_files);
mat_files = mat_files(idx);

v_fac = {'0.1', '0.03', '1.0', '3.0', '10.0'};
colors = parula(length(v_fac));

L = nan(size(v_fac));

figure
for j=1:length(v_fac)
    idx = cellfun(@(x) contains(x,['v' v_fac{j}]),mat_files);
    mat_filesi = mat_files(idx);


    
    for i=1:length(mat_filesi)
        % Plot results
        load(['../multiphase_test/' mat_filesi{i}]);
        
        A = out{1}; zvec = out{2}; pvec = out{3}; ugvec = out{4}; umvec = out{5};
        phivec = out{6}; rhogvec = out{7}; chidvec = out{8}; Qmvec = out{9}; Qgvec = out{10}; failure = out{11};
        
        lithvec = -(.25/(1-.25))*zvec*2700*9.8;
        
        shearstress = nan(size(umvec));
        ss_bbl = 8*A.mu(phivec,pvec).*umvec./A.r;
        ss_frag = A.f0.*rhogvec.*ugvec.^2/A.r;

        shearstress(zvec<A.fragdepth) = ss_bbl(zvec<A.fragdepth);
        shearstress(zvec>=A.fragdepth) = ss_frag(zvec>=A.fragdepth);
        
        L(j) = semilogx(shearstress./lithvec,zvec,'Color',colors(j,:)); hold on;
        %semilogx(pvec,zvec); hold on;
        ylabel('depth');
        xlabel('pressure');
    end
end

legend(L,v_fac)


%semilogx(pvec,zvec,'-k','LineWidth',3)
xlim([1e5, 5e8]);


%% Get one example to save as CSV
load(['../multiphase_test/v0.03T1100.0r50.0d6000.0.mat']);
A = out{1}; zvec = out{2}; pvec = out{3}; ugvec = out{4}; umvec = out{5};
phivec = out{6}; rhogvec = out{7}; chidvec = out{8}; Qmvec = out{9}; Qgvec = out{10}; failure = out{11};
        
        lithvec = -(.25/(1-.25))*zvec*2700*9.8;
        
        shearstress = nan(size(umvec));
        ss_bbl = 8*A.mu(phivec,pvec).*umvec./A.r;
        ss_frag = A.f0.*rhogvec.*ugvec.^2/A.r;

        shearstress(zvec<A.fragdepth) = ss_bbl(zvec<A.fragdepth);
        shearstress(zvec>=A.fragdepth) = ss_frag(zvec>=A.fragdepth);
        
        L(j) = semilogx(shearstress./lithvec,zvec,'Color','r'); hold on;
        %semilogx(pvec,zvec); hold on;
        ylabel('depth');
        xlabel('\tau');
        
        pfile = [zvec, pvec];
        shearfile = [zvec, shearstress];
        csvwrite('p_v0.03.csv',pfile);
        csvwrite('tau_v0.03.csv',shearfile);

