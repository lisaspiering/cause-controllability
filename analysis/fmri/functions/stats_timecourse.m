function [all_peaks, all_peaks_time, all_stats] = stats_timecourse(allBeta,modelMap,cfg, ev, reg)
% this function does the ttest against 0 for the current model
% we do leave-one-out procedure
% output is saved in an excel sheet that has separate sub-sheets for each
% regressor this function is run for
% Lisa Spiering, 2025

% some general parameters
rois = cfg.roi; % just copy over rois to test -> can also specify specific ones
rois_names = cfg.roi_names;

% Get data
currModel = modelMap(ev);
reg       = find(strcmp(currModel, reg));
numsubs   = 1:numel(cfg.subs2run); % get subject list as indices

% we want to run stats on timewindow 4-10sec
nsamples         = round(((cfg.pre_win+cfg.post_win)./cfg.TR)*cfg.upsample);
time             = -cfg.pre_win:(cfg.pre_win+cfg.post_win)/nsamples:cfg.post_win;
% find indices for test window
start_ind = find(time >= cfg.test_window(1),1,'first'); % 90 => time(90) = 4.0520s
end_ind   = find(time <= cfg.test_window(2),1,'last'); % 177 => time(177) = 9.9680s

peaks = []; peaks_time = []; % set up and clear content of these variables
for r = 1:numel(rois)
    for s = 1:length(numsubs) % iterate through participants
        numsubs2avg = numsubs; numsubs2avg(s) = []; % drop this participant from participants (i.e. leave-one-out)
        window  = mean(squeeze(allBeta.(sprintf('%s',rois{r})).(sprintf('%s',ev))(reg,:,numsubs2avg)),2); % avg across all timepoints for all remaining participants except the current one
        % find max of abs activation within the window
        [~,i] = max(abs(window(start_ind:end_ind)));
        i = i + (start_ind-1); % LEAVE! (adds peak to the overall time window of time courses (so that we have index not only within the timewindow 4-10sec
        peaks_time(s,r) = time(i);
        peaks(s,r) = allBeta.(sprintf('%s',rois{r})).(sprintf('%s',ev))(reg,i,s); % now takes peak at the timepoint for this sub
        clear window
    end % end participants
end % end rois

% do stats
T_pval = []; T_stats = []; T_DF = [];MU = []; SD = []; ROI = []; ROI_name = []; REG = []; 
for r = 1:numel(rois)
    % do onesample t test against 0
    [~,p,~,stats]  = ttest(peaks(:,r)); % test against 0 by default
    T_pval         = [T_pval; p];
    T_stats        = [T_stats; stats.tstat];
    T_DF           = [T_DF; stats.df];
    % general descriptive
    MU             = [MU; mean(peaks(:,r))];
    SD             = [SD; stats.sd];
    ROI            = [ROI; {rois{r}}];
    ROI_name       = [ROI_name; {rois_names{r}}];
    REG            = [REG; currModel(reg)];
end

% aggregate data into 3 tables
all_stats = table(REG, ROI_name,ROI,MU,SD,T_stats,T_DF,T_pval) %, W_zval, W_pval
all_peaks = array2table(peaks,'VariableNames',rois); all_peaks.ID = cfg.subs2run'; all_peaks.REG = repmat(currModel(reg),length(numsubs),1);% table with participant peak values per roi
all_peaks_time = array2table(peaks_time,'VariableNames',rois); all_peaks_time.ID = cfg.subs2run'; % table with participant peak times per roi

% save stats to excel sheet
suffix = ['_',num2str(cfg.test_window(1)),'to',num2str(cfg.test_window(2))]; % suffix for tested time window
writetable(all_stats,fullfile(cfg.roi_dir,['roi02c_stats_glm',cfg.glmix,'_',ev,suffix,'.xlsx']), 'Sheet', reg); %write table of this regressor into a separate sheet of the excel file
writetable(all_peaks,fullfile(cfg.roi_dir,['roi02c_stats_glm',cfg.glmix,'_',ev,suffix,'_subpeaks.xlsx']), 'Sheet', reg);
end