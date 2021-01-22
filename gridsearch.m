% gridsearch.m
% Used to explore the solution space for the muphem model

%% CONDUCT GRID SEARCH

% Get list of json input files
json_dir = '../multiphase_test2';
jsons = dir(json_dir);
jsons = {jsons.name};
idx = cellfun(@(x) contains(x,'.json'),jsons);
jsons = jsons(idx);

A = Amodels.initA_paper_mod;

for i=1:length(jsons)
    % Run muphem code
    input_json = ['../multiphase_test2/' jsons{i}];
    params = jsondecode(fileread(input_json));
    try
        out = muphem('multiflow2',0,A,params);

        % Save output
        output_name = [input_json(1:end-5) '.mat'];
        save(output_name, 'out');
        disp('Made a save!')
    catch ERROR
        disp(ERROR)
        
    end

end

return

%% xc
mat_dir = './multiphase_test';
mat_files = dir(mat_dir);
mat_files = {mat_files.name};
idx = cellfun(@(x) contains(x,'.mat'),mat_files);
mat_files = mat_files(idx);

xc_fac = { '0.0',  '0.2', '0.4', '0.6'};
colors = parula(length(xc_fac));

L = nan(size(xc_fac));

figure
for j=1:length(xc_fac)
    idx1 = cellfun(@(x) contains(x,['xc' xc_fac{j}]),mat_files);
    idx2 = cellfun(@(x) contains(x,['d5500']),mat_files);
    idx3 = cellfun(@(x) contains(x,['T1300']),mat_files);
    idx = idx1 & idx2 & idx3;
    mat_filesi = mat_files(idx);


    
    for i=1:length(mat_filesi)
        % Plot results
        load(['./multiphase_test/' mat_filesi{i}]);
        
        A = out{1}; zvec = out{2}; pvec = out{3}; ugvec = out{4}; umvec = out{5};
        phivec = out{6}; rhogvec = out{7}; chidvec = out{8}; Qmvec = out{9}; Qgvec = out{10}; failure = out{11};
        
        %lithvec = -(.25/(1-.25))*zvec*2700*9.8;
        lithvec = -zvec*2700*9.8;
        
        shearstress = nan(size(umvec));
        ss_bbl = 8*A.mu(phivec,pvec).*umvec./A.r;
        ss_frag = A.f0.*rhogvec.*ugvec.^2/A.r;

        shearstress(zvec<A.fragdepth) = ss_bbl(zvec<A.fragdepth);
        shearstress(zvec>=A.fragdepth) = ss_frag(zvec>=A.fragdepth);
        
        %L(j) = semilogx(shearstress./lithvec,zvec,'Color',colors(j,:)); hold on;
        L(j) = semilogx(pvec./lithvec,zvec,'Color',colors(j,:)); hold on;
        %semilogx(pvec,zvec); hold on;
        xlim([1e-2, 1e2]);
        ylabel('depth');
        xlabel('pressure/lithostatic');
        title('\chi_c')
    end
end

legend(L,xc_fac)

return

%% radius
mat_dir = './multiphase_test';
mat_files = dir(mat_dir);
mat_files = {mat_files.name};
idx = cellfun(@(x) contains(x,'.mat'),mat_files);
mat_files = mat_files(idx);

xc_fac = { '1000',  '1100', '1200', '1300'};
colors = parula(length(xc_fac));

L = nan(size(xc_fac));

figure
for j=1:length(xc_fac)
    idx1 = cellfun(@(x) contains(x,['T' xc_fac{j} '.0r']),mat_files);
    idx2 = cellfun(@(x) contains(x,['d5500']),mat_files);
    idx3 = cellfun(@(x) contains(x,['xc0.2']),mat_files);
    idx = idx1 & idx2 & idx3;
    
    mat_filesi = mat_files(idx);


    
    for i=1:length(mat_filesi)
        % Plot results
        load(['./multiphase_test/' mat_filesi{i}]);
        
        A = out{1}; zvec = out{2}; pvec = out{3}; ugvec = out{4}; umvec = out{5};
        phivec = out{6}; rhogvec = out{7}; chidvec = out{8}; Qmvec = out{9}; Qgvec = out{10}; failure = out{11};
        
        %lithvec = -(.25/(1-.25))*zvec*2700*9.8;
        lithvec = -zvec*2700*9.8;
        
        shearstress = nan(size(umvec));
        ss_bbl = 8*A.mu(phivec,pvec).*umvec./A.r;
        ss_frag = A.f0.*rhogvec.*ugvec.^2/A.r;

        shearstress(zvec<A.fragdepth) = ss_bbl(zvec<A.fragdepth);
        shearstress(zvec>=A.fragdepth) = ss_frag(zvec>=A.fragdepth);
        
        %L(j) = semilogx(shearstress./lithvec,zvec,'Color',colors(j,:)); hold on;
        L(j) = semilogx(pvec./lithvec,zvec,'Color',colors(j,:)); hold on;
        drawnow
        %semilogx(pvec,zvec); hold on;
        ylabel('depth');
        xlabel('pressure/lithostatic');
        title('temperature')
    end
end

legend(L,xc_fac)

%semilogx(pvec,zvec,'-k','LineWidth',3)
xlim([1e-2, 1e2]);


%% Get one example to save as CSV
load(['./multiphase_test/xc0.2T1300.0r10.0d5500.0.mat']);
A = out{1}; zvec = out{2}; pvec = out{3}; ugvec = out{4}; umvec = out{5};
phivec = out{6}; rhogvec = out{7}; chidvec = out{8}; Qmvec = out{9}; Qgvec = out{10}; failure = out{11};
        
        lithvec = -zvec*2700*9.8;
        
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
        csvwrite('p_T1300.csv',pfile);
        csvwrite('tau_T1300.csv',shearfile);

