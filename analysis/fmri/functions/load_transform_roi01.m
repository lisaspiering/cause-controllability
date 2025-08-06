function dataB = load_transform_roi01(cfg)

% load behavioural data (is output from cc_roi01 R script)
tmp = load(fullfile(cfg.roi_dir,['roi01_behdata_glm' cfg.glmix '.mat']));

% transform data to table manually bc somehow readtable just doesn't work anymore
dataB = table();
for i = 1:size(tmp.namesB,1)
    % sometimes the columns need to be transposed, sometimes not
    [r, c] = size(tmp.dataB.(tmp.namesB{i})); % check size
    if r == 1
        dataB.(tmp.namesB{i}) = double(tmp.dataB.(tmp.namesB{i})');
    elseif c == 1
        dataB.(tmp.namesB{i}) = double(tmp.dataB.(tmp.namesB{i}));
    else
        disp(['problem: check column dimensions for variable: ' tmp.namesB{i}])
    end
end
dataB.eoi = cellstr(tmp.eois); % manually add column that contains strings, somehow for this design 03 we need to convert manually from character to cell array
end