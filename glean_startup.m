here = fileparts(mfilename('fullpath'));

% Add OHBA shared libraries
if ~exist(fullfile(here,'..','osl-external2'))
    fprintf(2,'Could not find ''%s''\n',fullfile(osldir,'osl-external'));
    error('osl-external is missing. Clone https://github.com/OHBA-analysis/osl-external into the same directory as osl2');
end

addpath(fullfile(here,'..','osl-external'));
osl_external_startup
