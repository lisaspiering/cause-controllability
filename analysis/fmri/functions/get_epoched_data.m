function [dataBupd, cfg] = get_epoched_data(dataB,cfg)
    % this functions takes the behavioural data that was created in roi01 R
    % script. For each participant, it epochs the time course data of each ROI
    % around each event of interest. Trials where an event of interest is 
    % too close to the end of the scan (this is considering the duration of 
    % the epoch), are removed. 
    %  N Khalighinejad Adapted 07/10/19; H Trier adapted
    % L Spiering adapted 21/11/22
    % Input:
    %   - dataB: the behavioural data of roi01 R script
    %   - 
    % Outputs:
    %   - timecourse data for each participant, each ROI, each event of interest
    %   - dataBupd: dataB table where trials are removed where the onset of an event of
    %   interest + the post-onset duration of the epoch exceeds the scan
    
    dataBupd = dataB; dataBupd(:,:) = []; % new empty table where we store all trials that didn't have an onset+post_win exceeding the scan duration
    num_omitted_rows = 0;

    for r = 1:numel(cfg.roi)

        trial_data = [];

        for s = cfg.subs2run

            s    = s{1};
            sstr = ['s',num2str(s,'%02.f')]; % put together string format for folder

            %% get trial data

            % extracted ROI timecourse file
            cfg.data_path = fullfile(cfg.timecourse_dir,sstr,'masks_tc',['tc_' cfg.roi{r} '.txt']);

            fclose all;

            %% read in fMRI data

            % read in timecourse
            tc = dlmread(cfg.data_path)';

            %% upsample

            % number of samples in epoch
            nsamples = round(((cfg.pre_win+cfg.post_win)./cfg.TR)*cfg.upsample);

            % normalise timecourse
            tc       = zscore(tc);

            % provide polynomial form of the cubic spline interpolant to the data
            tc_x    = 1:length(tc);
            tc_yint = 1:(1/cfg.upsample):length(tc); % query points
            tc_int  = spline(tc_x,tc,tc_yint);

            %% Loop through events of interest
            for ev = 1:numel(cfg.eoi)

                e = cfg.eoi{ev};
                % get the event onset time from behavioural data
                evs_data  = dataB(strcmp(dataB.eoi,e) & dataB.ID==s,:); % filter by event of interest and current subject

                % get the onsets - depends on event of interest
                ev_onsets = evs_data.onset; % put together onset column e.g. onset_O for event O_p1

                %% epoch data
                onset = ev_onsets-cfg.pre_win-(cfg.TR/2);

                if onset(1) <= 0
                    error ('onset error')
                end

                % get the upsampled position of each onset
                onset = round(((onset)./cfg.TR)*cfg.upsample);

                % get indices into timecourse for each trial
                epoch_idx = repmat(onset,1,nsamples)+repmat(0:nsamples-1,size(onset,1),1);

                % Remove any trial that is close enough to the end that its
                % indices exceed the length of the timecourse
                omit_rows = [];
                for trl = 1:size(epoch_idx,1)
                    exceeds_tc = find(epoch_idx(trl,:)>size(tc_int,2));
                    if ~isempty(exceeds_tc)
                        if r==1 % display info of omitted trial only for first roi, don't need it again afterwards
                            disp(strcat(sstr," - trial ",num2str(trl)," out of ",num2str(size(epoch_idx,1))," trials exceeds timecourse / event: ",e))
                        end
                        omit_rows = [omit_rows, trl];
                    end
                end
                if r==1 % store number of omitted rows only for first roi
                    num_omitted_rows = num_omitted_rows + length(omit_rows);
                end
                epoch_idx(omit_rows,:) = [];

                % update the removed events also in behaviour (and then save
                % updated behavioural file in a new file (but only if we're in the
                % first ROI, because we only need to do this once)
                evs_data(omit_rows, :) = [];
                if r == 1 % only for first ROI that doesn't differ for ROIs
                    dataBupd = [dataBupd;evs_data];
                end

                % Index timecourse
                trial_data = tc_int(epoch_idx);

                % write epoched data
                trial_data = tc_int(epoch_idx);
                if  cfg.textout
                    % save data in subfolder for roi
                    idir = fullfile(cfg.timecourse_dir,sstr,['tc_glmix' cfg.glmix],cfg.roi{r});
                    % create subfolder if doesn't exist yet
                    if ~exist(idir, 'dir') 
                        mkdir(idir)
                    end
                    filename = fullfile(idir,['tc_' cfg.roi{r},'_',e,'_epoched']);
                    save(filename,'trial_data');
                end

                clear trial_data epoch_idx onset idir
            end

            clear tc tc_int tc_yint tc_x
        end

    end

    % save updated behavioural data table
    % this is the behavioural data without the events that exceeded the time
    if cfg.textout
        fileoutname = ['roi02a_epoched_data_glm' cfg.glmix '.mat'];
        save(fullfile(cfg.roi_dir,fileoutname),'dataBupd','cfg')
    end

end