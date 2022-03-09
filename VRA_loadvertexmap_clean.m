function Tout = VRA_loadvertexmap_clean(vectordata)

% Load and clean data
% load raw data, and clean data by replacing 0, NaN and outliers
% Input
% ---- vectordata: vectordata to be preprocessed
% Output
% ---- Tout:     the Matlab table containing the data after cleaning

T = vectordata;


%% Replace cortical thickness 0 with NaN
fprintf('Replacing all cortical thickness 0s with NaN\n');
T(T==0)  = NaN; % the cortical thickness data are in rows


%% Replace cortical thickness outliers (beypnd column mean +/- N*SD) with NaN
Ncoff = 3; % outliers should beyond mean +/- Ncoff*SD, default is 3, but lead to a lot of outliers in this data
fprintf('Replacing cortical thickness outliers (beyond column mean +/- %d*SD) with NaN\n',Ncoff);
for j = 1:size(T,1) % per row of CT
    T1 = T(j,:);
    Tmean = nanmean(T1);
    Tstd  = nanstd(T1);
    % look for upper outliers
    idx1 = []; idx2 = []; idx = [];
    idx1 = find(T1 > Tmean + Ncoff*Tstd); % index of the upper outliers
    idx2 = find(T1 < Tmean - Ncoff*Tstd); % index of the lower outliers
    idx  = [idx1,idx2];
    
    if idx % if there is at least an outlier
        % display the SubjID and cortical area name of the outlier
        fprintf('--Outliers detected in vertex %d, numsubs = %d', j, length(idx));
        %for i = 1:length(idx), fprintf('%d(%1.3f<>%1.3f+-%d*%1.3f),',T.output.subject_names{idx(i)},j,Tmean,Ncoff,Tstd); end
        
        fprintf('\n');
        % replace outliers
        T(j,idx) = NaN;
    end
end
fprintf('Overal values: Min=%1.3f, Mean=%1.3f, Max=%1.3f\n', min(nanmin(T)), mean(nanmean(T)), max(nanmax(T)));

%% Replace cortical thickness NaN with Column Mean per group
fprintf('Replacing cortical thickness NaN with Column Mean per group\n');
T1 = T;
for j = 1:size(T,1) % per column
    col_mean = nanmean(T1(j,:));
    if isnan(col_mean)
        col_mean = 0;
    end
    T1(j,isnan(T1(j,:))) = col_mean;
end
T = T1;


%% Output
Tout = T;
fprintf('Number of subjects: %d\n', size(Tout,2));
fprintf('==============\n==============\n==============\nLoad and Clean Completed!!!\n\n\n');

end