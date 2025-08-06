function plot_timecourse(allBeta,modelMap,cfg)
% this function plots the timecourses for the previously epoched + glm-ed
% data.
% Lisa Spiering, 2025

% plotting variables
nsamples = round(((cfg.pre_win+cfg.post_win)./cfg.TR)*cfg.upsample);
time     = -cfg.pre_win:(cfg.pre_win+cfg.post_win)/nsamples:cfg.post_win;

% color for beta effect
plt_colour = {{'red',[0.8,0,0]},{'blue',[0,0,0.8]}};

% loop through events of interest
for ev = cfg.eoi
    ev = ev{1};
    currModel = modelMap(ev);
    for contrast = 1:length(currModel)
        z = 1;
        if numel(cfg.roi) <= 5
            figure('Position',[500 500 1000 300],'color','w'); % position: [left bottom width height]
            sp_rows = 1; sp_cols = numel(cfg.roi); % 1 row for plots if < 4 rois, as many cols as rois
        else
            figure('Position',[500 500 1000 500],'color','w'); 
            sp_rows = 2; sp_cols = 4;%round(numel(cfg.roi)/3);
        end
        
        for r = 1:numel(cfg.roi)
            % generate subplot and general settings
            subplot(sp_rows, sp_cols,z);
            hold on
            
            % compute betas for this contrast, roi, ev
            beta = squeeze(allBeta.(sprintf('%s',cfg.roi{r})).(sprintf('%s',ev))(contrast,:,:));
            
            % plot beta with sem shading
            semshade_LS(beta',0.15,plt_colour{1}{2},linspace(-cfg.pre_win,cfg.post_win,nsamples)');

            % add another beta in a different color if specified
            if contrast==cfg.regs_in_one(1)
                beta = squeeze(allBeta.(sprintf('%s',cfg.roi{r})).(sprintf('%s',ev))(cfg.regs_in_one(2),:,:)); % compute beta for 2nd contrast to plot
                semshade_LS(beta',0.15,plt_colour{2}{2},linspace(-cfg.pre_win,cfg.post_win,nsamples)'); % plot in 2nd color
            end
            
            % set ylims higher for constant than for others
            if strcmp(currModel{contrast},'constant')
                ylim([-0.6 0.6]);
            else
                ylim([-0.25 0.25]);
            end
            
            xlim([(-cfg.pre_win) cfg.post_win]);
            set(gca,'fontsize',20);
            yL = get(gca,'YLim'); line([0 0],yL,'Color',[0.8 0.8 0.8]);
            xL = get(gca,'XLim'); line(xL,[0 0],'Color',[0.8 0.8 0.8]);
            % Label plot
            title(strrep(cfg.roi_names{r},' [','\fontsize{18} [')); xlabel(gca, 'Time (s)'); %ylabel(gca, 'Effect on BOLD signal'); %
            
            z = z+1;
        end
        % title
        if contrast==cfg.regs_in_one(1)
            sgtitle([plt_colour{1}{1} ' ' currModel{contrast} ' and ' plt_colour{2}{1} ' ' currModel{cfg.regs_in_one(2)} ' during ' ev], 'Interpreter', 'none','FontWeight','bold','FontSize',25)
        else 
            sgtitle([currModel{contrast} ' during ' ev], 'Interpreter', 'none','FontWeight','bold','FontSize',25)
        end
    end
end
end