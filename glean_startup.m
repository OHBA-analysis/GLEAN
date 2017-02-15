here = fileparts(mfilename('fullpath'));

% Add OHBA shared libraries
if ~exist(fullfile(here,'..','ohba-external'))
    fprintf(2,'Could not find ''%s''\n',fullfile(here,'..','ohba-external'));
    error('ohba-external is missing. Clone https://github.com/OHBA-analysis/ohba-external into the same directory as GLEAN');
end

addpath(fullfile(fileparts(here),'ohba-external'));
ohba_external_startup
