function [allBeta, modelMap] = run_timecourse_glm(dataBupd,modelMap,cfg)
% this script runs the specified GLM
% 21/11/22 L Spiering adapted from H Trier's version of N Khalighinejad script

%% Filter based on action type & set up model
for ev = cfg.eoi
    
    ev       = ev{1};
    tmp      = split(ev,'_p');
    
    % Get data for this event - remove task phases we don't need
    dataBev   = dataBupd(strcmp(dataBupd.eoi,ev),:);
    
    % loop through ROI
    for r = 1:numel(cfg.roi)
        z = 0; % separate counter again for each person
        for s = cfg.subs2run
            s = s{1};
            sstr = ['s',num2str(s,'%02.f')]; % put together string format for folder
            z = z+1;  % counter for each person
            
            dataBs = dataBev(dataBev.ID==s, :); % get behavioural data of this participant
            
            
            % load time-series data (variable trial_data)
            tc_dir = fullfile(cfg.timecourse_dir,sstr,['tc_glmix' cfg.glmix],cfg.roi{r});
            load(fullfile(tc_dir,['tc_',cfg.roi{r},'_',ev,'_epoched']));
            
            %% run GLM
            % regressor design matrix (don't forget constant!)
            dmat = table2array(dataBs(:,modelMap(ev)));
            
            % do we need to normalise the regressors here or did we already
            % do that in R script before?
            mustscale = dataBs.mustscale(1);
            
            % contrast design matrix
            contrasts = diag(ones(size(dmat,2),1));
            
            % normalise data - don't change position of this! / it's important that constant goes last for this! (so that this is not normalised)
            % update 07.03.2023 LS: we custom-normalise
            % regressors and check whether we normalise here again with a
            % helper variable in dmat
            if mustscale
                dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1)));
            end
            
            % beta X time output
            if size(trial_data,1) == size(dmat,1)
                betaOut = ols(trial_data,dmat,contrasts); % regression
                allBeta.(sprintf('%s',cfg.roi{r})).(sprintf('%s',ev))(:,:,z) = betaOut; % save beta in overall struct
                
                % write effect size / save individual beta file
                if  cfg.textout
                    filename = fullfile(tc_dir,[cfg.roi{r},'_',ev,'_effectSize']);
                    save(filename,'betaOut');
                end
            else
                disp(strcat("Alert! Trial data and dmat not equal sizes for event ",ev));
            end
            fclose all;
            clear dmat contrasts trial_data betaOut dataBs
        end
    end
    clear dataBev
end

% save allBeta
if cfg.textout
    save([cfg.roi_dir 'roi02b_tc_glm' cfg.glmix],'allBeta','modelMap');
    save([cfg.roi_dir 'roi02b_tc_glm' cfg.glmix '_allBeta'],'allBeta'); % for R 
end
end