here = fileparts(mfilename('fullpath'));

% Add OHBA shared libraries
if ~exist(fullfile(here,'..','osl-external'))
    fprintf(2,'Could not find ''%s''\n',fullfile(here,'..','osl-external'));
    error('osl-external is missing. Clone https://github.com/OHBA-analysis/osl-external into the same directory as GLEAN');
end

addpath(fullfile(fileparts(here),'osl-external'));
osl_external_startup
