function pathname = glean_path
% Returns the path to the GLEAN code.
%
% pathname = GLEAN_PATH
%
% Adam Baker 2015

[pathname,~,~] = fileparts(which('glean_run'));

end
