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

idx = cellfun(@(x) contains(x,'v0.03'),mat_files);
mat_files = mat_files(idx);


figure
for i=1:length(mat_files)
    % Plot results
    load(['../multiphase_test/' mat_files{i}]);
    
    A = out{1}; zvec = out{2}; pvec = out{3}; ugvec = out{4}; umvec = out{5};
    phivec = out{6}; rhogvec = out{7}; chidvec = out{8}; Qmvec = out{9}; Qgvec = out{10}; failure = out{11};
   
    semilogx(pvec,zvec); hold on;
    
end
