function GLEAN = data(GLEAN)
% Set up the directory structure for a new or existing GLEAN analysis.
%
% GLEAN = glean.data(GLEAN)
%
% Adam Baker 2015


% Check data exists and is in the right format
if ~isfield(GLEAN,'data') || ~all(cellfun(@ischar,GLEAN.data))
    error('Must specify GLEAN.data as a [sessions x 1] cell array"')
end
   
% Loop through each module and set up the directory and file structure
for module = {'envelope','subspace','model'}
    GLEAN = setup_files(GLEAN,module);
end

end


function GLEAN = setup_files(GLEAN,module)
% Set up the list of files for each stage
[~,sessionNames] = cellfun(@fileparts,GLEAN.data,'UniformOutput',0);

GLEAN.(char(module)).data = fullfile(GLEAN.(char(module)).dir,'data',strcat(sessionNames,'.mat'));
if ~isdir(fullfile(GLEAN.(char(module)).dir,'data'))
    mkdir(fullfile(GLEAN.(char(module)).dir,'data'));
end

if strcmp(module,'model')
    GLEAN.model.model = fullfile(GLEAN.model.dir,'model.mat');
end

end




