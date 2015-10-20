function GLEAN = glean_check(GLEAN)
% Checks the settings and set up directory and data structures in GLEAN.
%
% GLEAN = GLEAN_CHECK(GLEAN)
%
% Adam Baker 2015

GLEAN = glean_settings(GLEAN);
GLEAN = glean_directories(GLEAN);
GLEAN = glean_data(GLEAN);