clear;
clc;

%% Inputs
rundir_path = '/home/viraj/Desktop/lab_morey/ncanda/exp_glm/baselinenobinge_hcombat_drkclss/run1';
vertexdata_filepath = fullfile('/home/viraj/Desktop/lab_morey/ncanda/exp_glm/baselinenobinge_hcombat_drkclss', ...
    'hComBat_thickness_fsaverage5_fwhm20_d20484_n2809.mat');

EOI = 'all';
num_comps = 7;

%
fprintf('NMF run: \n\t%s\n',rundir_path);
fprintf('No. of components: %d\n', num_comps);

%% setup paths
components_filepath = fullfile(rundir_path, 'Output', ['Result_',EOI,'_',sprintf('%03d',num_comps),'_components.mat']);

output_filepath = fullfile('/home/viraj/Desktop/lab_morey/ncanda/exp_glm/baselinenobinge_hcombat_drkclss/glm_all/', ...
    ['NMF','_k', sprintf('%02d',num_comps), '_', 'd20484_n2809', '.mat']);

output_xclean_filepath = fullfile('/home/viraj/Desktop/lab_morey/ncanda/exp_glm/baselinenobinge_hcombat_drkclss/glm_all/', ...
    ['XClean_', 'd20484_n2809', '.mat']);

%% Load NMF components
fprintf('Loading components from : \n\t%s\n',components_filepath);
nmf_data = load(components_filepath);
W = nmf_data.W;

%% Preprocessing the raw vertex data
fprintf('Processing data for NMF ...');
fprintf('Loading data from : \n\t%s\n',vertexdata_filepath);
data = load(vertexdata_filepath);

% load relevant subjects
subjects = data.output.subject_names;

% Functions
% Clean data

disp('Preprocessing using vertex pipeline');
X_processed = VRA_loadvertexmap_clean(data.output.vectordata);  

%% Obtaining coefficients of NMF components for subjects

H = W'*X_processed;


num_subs = size(H, 2);
subjectids = subjects;


fprintf('NMF subjects: %d\n',num_subs);


%% Removing overlaps in the components
% cleaning W using the class map

K_clusters = clusteringAssignment(W);

W_clean = zeros(size(W));


for i=1:num_comps
    W_clean(K_clusters == i, i) = W(K_clusters==i, i);
end

%% Get H in physical units using Aris' function

H_physical = AUtoPhysicalUnits(X_processed, W_clean);

%% Save H_physical, W_clean, K_clusters, X_processed, H, W

save(output_xclean_filepath, 'X_processed', 'subjectids');
save(output_filepath, 'H_physical', 'W_clean', 'K_clusters', 'H', 'W', 'subjectids');
