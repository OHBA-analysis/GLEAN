function GLEAN = load(glean_file)
% Load an existing GLEAN analysis
%
% GLEAN = glean.load(glean_file)
%
% REQUIRED INPUTS:
%   GLEAN     - A .mat file containing an existing GLEAN analysis
%
% Adam Baker 2015

[p,f,e] = fileparts(glean_file);

if isempty(e)
    e = '.mat';
end

if isempty(p)
    p = pwd;
end
    
if ~isequal(e,'.mat')
    error('Not a valid glean file');
end

glean_file = fullfile(p,[f,e]);


try 
    tmp = load(glean_file);
    GLEAN = tmp.GLEAN;
catch
    error('Not a valid glean file');
end

GLEAN.name = glean_file;

GLEAN = glean.check(GLEAN);

save(GLEAN.name,'GLEAN')

end