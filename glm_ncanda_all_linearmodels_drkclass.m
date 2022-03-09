clear;
clc;

%% Inputs
rundir_path = './';
demogs_path = './ncanda_all_sites_out_drkclss_20200617.xlsx';

EOI = 'all';
num_comps = 10;
analysis_name = 'maineffects';
bonferroni = true;
pval_thresh = 0.05;

%
fprintf('NMF run: \n\t%s\n',rundir_path);
fprintf('No. of components: %d\n', num_comps);
fprintf('GLM name: %s\n', analysis_name);
%% setup paths
filepath = fullfile(rundir_path, ['Result_',EOI,'_',sprintf('%03d',num_comps),'_components.mat']);
data_filepath = fullfile(rundir_path,'hComBat_thickness_fsaverage5_fwhm20_d20484_n2809.mat');

%% Computed required Statistics before subset selection
% Load demographics
fprintf('Loading demographics from : \n\t%s\n', demogs_path);
demogs = get_ncanda_demogs_drkclass(demogs_path);

fprintf('Computing required metrics:\n');
% log(ICV)
fprintf('\tlog(ICV)\n');
demogs.logicv = log(demogs.ICV);

% age_d
subject_ids = unique(demogs.subject);
sub_mean_age = zeros(size(demogs.mri_t1_age));
for i = 1:size(subject_ids,1)
    subj_idx = demogs.subject == subject_ids(i,:);
    sub_mean_age(subj_idx) = mean(demogs.mri_t1_age(subj_idx));
end

% binged_YN
num_binged = 0;
num_nonbinged = 0;
binged_yn = zeros(size(demogs.mri_t1_age));
for i = 1:size(subject_ids,1)
    subj_idx = demogs.subject == subject_ids(i,:);
    cddrsum = sum(demogs.cddr_past_year_binge(subj_idx));
    if cddrsum > 0
        binged_yn(subj_idx) = 1;
        num_binged = num_binged+1;
    else
        num_nonbinged = num_nonbinged+1;
    end
end
fprintf('\tnum_binged: %d\n',num_binged);
fprintf('\tnum_nonbinged: %d\n',num_nonbinged);
demogs.binged_yn = binged_yn;

fprintf('\tage_d\n');
demogs.age_d = demogs.mri_t1_age - sub_mean_age;
fprintf('\tage_m\n');
demogs.age_m = sub_mean_age - mean(sub_mean_age);


%% Load NMF components
fprintf('Loading components from : \n\t%s\n',filepath);
nmf_data = load(filepath);
W = nmf_data.W;
H = nmf_data.H;


%% Preprocessing the raw vertex data
fprintf('Processing data for NMF ...');
fprintf('Loading data from : \n\t%s\n',data_filepath);
data = load(data_filepath);

% load relevant subjects
subjects = data.output.subject_names;

% Functions
% Clean data

disp('Preprocessing using vertex pipeline');
X = VRA_loadvertexmap_clean(data.output.vectordata);  

%% Obtaining coefficients for subjects
H = W'*X;


num_subs = size(H, 2);
nmf_subjectids = subjects;


fprintf('NMF subjects: %d\n',num_subs);
%% Create a data frame
% syncing ids in demogs to nmf_subjects
fprintf('Selecting demographics for NMF subjects\n');
[~,demog_idxs] = ismember(nmf_subjectids,demogs.processed_id);
data = demogs(demog_idxs,:);

%%

% removing outliers or invalid rows.
analysis_idxs = ~(data.ICV == 0);
data_analysis = data(analysis_idxs, :);
H_analysis = H(:, analysis_idxs);

if bonferroni
    pval_threshold = pval_thresh/num_comps;
    fprintf('Bonferroni correction applied\n');
else
    pval_threshold = pval_thresh;
end

fprintf('Pval threshold: %d\n', pval_threshold);
% 
% checkp_name = {'exceeds_bl_drinking_Y', ...
%     'cddr_past_year_binge', ...
%     'fh_alc_density', ...
%     'life_trauma_RP',...
%     'age_d:age_m', ...
%     'cddr_past_year_binge:age_d:age_m', ...
%     'cddr_past_year_binge:age_d', ...
% 	'cddr_past_year_binge:age_m', ...
% 	'cddr_past_year_binge:life_trauma_RP:age_d', ...
%     'life_trauma_RP:age_d:age_m'};
% 
checkp_name = {'exceeds_bl_drinking_Y', ...
    'cddr_past_year_binge', ...
    'fh_alc_density', ...
    'life_trauma_RP',...
    'drinking_class',...
    'mri_t1_age:drinking_class', ...
    'life_trauma_RP',...
    'exceeds_bl_drinking_Y:age_d',...
    'exceeds_bl_drinking_Y:age_m',...
    'drinking_class:age_d',...
    'drinking_class:age_m',...
    'drinking_class:age_d:age_m',...
    'life_trauma_RP:drinking_class:age_d',...
    'life_trauma_RP:drinking_class',...
    'life_trauma_RP:age_d'};

%LOGICV
% models = {'nmf_component ~ age_d + age_m + exceeds_bl_drinking + logicv + sex + ses + (1|participant_id) + (1|site)', ...
%     'nmf_component ~ age_d + age_m + cddr_past_year_binge + logicv + sex + race + ses + fh_alc_density + life_trauma_RP + (1|participant_id) + (1|site)', ...
%     'nmf_component ~ age_d*age_m*cddr_past_year_binge + logicv + sex + race + ses + fh_alc_density + life_trauma_RP + (1|participant_id) + (1|site)', ...
%     'nmf_component ~ age_d*age_m*life_trauma_RP +cddr_past_year_binge + logicv + sex + race + ses + fh_alc_density + (1|participant_id) + (1|site)', ...
%     'nmf_component ~ age_d*life_trauma_RP*cddr_past_year_binge + age_m + logicv + sex + race + ses + fh_alc_density + (1|participant_id) + (1|site)'};
% 
% No logicv
% models = {'nmf_component ~ age_d + age_m + exceeds_bl_drinking + sex + ses + (1|participant_id) + (1|site)', ...
%     'nmf_component ~ age_d + age_m + cddr_past_year_binge + sex + race + ses + fh_alc_density + life_trauma_RP + (1|participant_id) + (1|site)', ...
%     'nmf_component ~ age_d*age_m*cddr_past_year_binge + sex + race + ses + fh_alc_density + life_trauma_RP + (1|participant_id) + (1|site)', ...
%     'nmf_component ~ age_d*age_m + cddr_past_year_binge + sex + race + ses + fh_alc_density + life_trauma_RP + (1|participant_id) + (1|site)', ...
%     'nmf_component ~ age_d*age_m*life_trauma_RP +cddr_past_year_binge + sex + race + ses + fh_alc_density + (1|participant_id) + (1|site)', ...
%     'nmf_component ~ age_d*life_trauma_RP*cddr_past_year_binge + age_m + sex + race + ses + fh_alc_density + (1|participant_id) + (1|site)'};

% mri_t1_age and no age_d and age_m
% models = {'nmf_component ~ mri_t1_age + exceeds_bl_drinking + sex + ses + (1|participant_id) + (1|site)', ...
%     'nmf_component ~ mri_t1_age + cddr_past_year_binge + sex + race + ses + fh_alc_density + life_trauma_RP + (1|participant_id) + (1|site)', ...
%     'nmf_component ~ mri_t1_age*cddr_past_year_binge + sex + race + ses + fh_alc_density + life_trauma_RP + (1|participant_id) + (1|site)', ...
%     'nmf_component ~ mri_t1_age*life_trauma_RP +cddr_past_year_binge + sex + race + ses + fh_alc_density + (1|participant_id) + (1|site)', ...
%     'nmf_component ~ mri_t1_age*life_trauma_RP*cddr_past_year_binge + sex + race + ses + fh_alc_density + (1|participant_id) + (1|site)'};

% models = {'nmf_component ~ mri_t1_age + exceeds_bl_drinking + sex + ses + (1|participant_id) + (1|site)', ...
%     'nmf_component ~ mri_t1_age + cddr_past_year_binge + sex + race + ses + fh_alc_density + life_trauma_RP + (1|participant_id) + (1|site)', ...
%     'nmf_component ~ mri_t1_age*cddr_past_year_binge + sex + race + ses + fh_alc_density + life_trauma_RP + (1|participant_id) + (1|site)'};

models = {'nmf_component ~ mri_t1_age + exceeds_bl_drinking + sex + ses + (1|participant_id)', ...
    'nmf_component ~ age_d + age_m + exceeds_bl_drinking + sex + ses + (1|participant_id)', ...
    'nmf_component ~ age_m + age_d*exceeds_bl_drinking + sex + ses + (1|participant_id)', ...
    'nmf_component ~ age_d + age_m*exceeds_bl_drinking + sex + ses + (1|participant_id)', ...
    'nmf_component ~ mri_t1_age + drinking_class + sex + race + ses + fh_alc_density + life_trauma_RP + (1|participant_id)', ...
    'nmf_component ~ age_d + age_m + drinking_class + sex + race + ses + fh_alc_density + life_trauma_RP + (1|participant_id)', ...
    'nmf_component ~ age_m*drinking_class + age_d + sex + race + ses + fh_alc_density + life_trauma_RP + (1|participant_id)', ...
    'nmf_component ~ age_d*drinking_class + age_m + sex + race + ses + fh_alc_density + life_trauma_RP + (1|participant_id)', ...
    'nmf_component ~ age_d*drinking_class*life_trauma_RP + age_m + sex + race + ses + fh_alc_density + (1|participant_id)', ...
    'nmf_component ~ age_d*age_m*drinking_class + sex + race + ses + fh_alc_density + life_trauma_RP + (1|participant_id)'};

for k = 1:size(models, 2)
    model = models{k};
    fprintf('Model: %d\n\t%s\n', k, model);

    %fprintf('Checking significance for : %s\n', checkp_name);

    fprintf('GLM on components : \n');
    pvalues = zeros(size(checkp_name, 2), size(H_analysis, 1));
    for i = 1:size(H_analysis, 1)
        comp = H_analysis(i,:);
        data_analysis.nmf_component = comp';
        comp_model = fitlme(data_analysis, model);
        if i == 1
            fprintf('Model coefficient names:\n');
            for m = 1:size(comp_model.CoefficientNames, 2)
                fprintf('\t%s\n',comp_model.CoefficientNames{m});
            end
        end
        fprintf('component: %d\n',i);
        % display results
        
        for j = 1:size(checkp_name,2)
            check_name = checkp_name{j};
            coeff_idx = strcmp(check_name, comp_model.CoefficientNames);
            if sum(coeff_idx) > 0
                p_value = comp_model.Coefficients.pValue(coeff_idx);
                beta = comp_model.Coefficients.Estimate(coeff_idx);

                fprintf('\tvariable: %s, p_value: %d, Beta: %d', check_name, p_value, beta); 

                if p_value < pval_threshold
                    fprintf(', Significant!!');
                end
                fprintf('\n');
                pvalues(j,i) = p_value;
            end
        end
    end
    
    fprintf('Component variable significance order:\n');
    for j = 1: size(pvalues, 1)
        if sum(pvalues(j, :)~=0) > 0
            fprintf('\t%s:\n', checkp_name{j})
            [sorted_pvals, pval_idxs] = sort(pvalues(j,:));
            for i = 1:size(pvalues,2)
                fprintf('\t%d',pval_idxs(i));
            end
            fprintf('\n')
            for i = 1:size(pvalues,2)
                fprintf('\t%d',sorted_pvals(i));
            end
            fprintf('\n');
        end
    end
    
    fprintf('\n');

end


%plotting for last model component vs mri_t1_age.
% decide which subjects did binging and which did not.
%binged_idx = data_analysis.cddr_past_year_binge > 0;
subject_ids_ana = unique(data_analysis.subject);
last_followup_idxs = zeros(size(subject_ids_ana));
for i = 1:size(subject_ids_ana,1)
    subj_idx = find(data_analysis.subject == subject_ids_ana(i,:));
    last_followup_idxs(i) = max(subj_idx);
end

selected_data = data_analysis(last_followup_idxs,:);
selected_h = H_analysis(:, last_followup_idxs);

binged_idx = selected_data.binged_yn > 0;

sig_components = [1,6,3,4,5,8];
xvals = data_analysis.mri_t1_age;%data_analysis.cddr_past_year_binge;%data_analysis.mri_t1_age.*data_analysis.cddr_past_year_binge;
for i = sig_components
    figure()
    plot(xvals(binged_idx), selected_h(i, binged_idx),'r.', 'markersize',10);
    hold on;
    lsline;
    plot(xvals(~binged_idx), selected_h(i, ~binged_idx),'b.', 'markersize',10);
    h = lsline;
    set(h, 'linewidth', 2);
    title(sprintf('component: %d',i));
    legend('binge drinking', '','controls','');
    xlabel('age');
    ylabel('coefficient');
end