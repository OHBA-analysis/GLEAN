function glean_timeseries(GLEAN)
% Runs the envelope stage of GLEAN.
%
% GLEAN_ENVELOPE(GLEAN)
%
% Adam Baker 2015


pretty_string('RUNNING TIMESERIES STAGE')

% Create temporary directory
[~,tmpdir] = fileparts(tempname);
tmpdir = fullfile(GLEAN.timeseries.dir,tmpdir);
mkdir(tmpdir);
c = onCleanup(@() system(['rm -r ' tmpdir]));

for session = 1:numel(GLEAN.data)

    % Check if envelope file exists and whether or not to overwrite
    file_exists = exist(GLEAN.timeseries.data{session},'file') == 2;
    overwrite   = GLEAN.timeseries.settings.overwrite == 1;

    subspace_method = char(intersect(fieldnames(GLEAN.subspace.settings),{'pca','parcellation','voxel'}));
    if strcmpi(subspace_method,'parcellation')
            msg = 'Subspace method is parcellation: skipping voxelwise timeseries computation\n';
            run_stage = false;
    elseif file_exists
        if overwrite
            msg = ['Overwriting existing timeseries file: \n' GLEAN.timeseries.data{session} '\n'];
            run_stage = true;
        else
            msg = ['Using existing envelope file: \n' GLEAN.timeseries.data{session} '\n'];
            run_stage = false;
        end
    else
        msg = ['Creating new envelope file: \n' GLEAN.timeseries.data{session} '\n'];
        run_stage = true;
    end
    fprintf(msg);


    if run_stage

        % Make a temporary filename to copy raw data to
        [~,tempdata] = fileparts(tempname);
        tempdata = fullfile(tmpdir,tempdata);

        % Copy data to temporary filename
        copymeeg(GLEAN.data{session},tempdata)

        % Compute envelopes:
        S               = [];
        S.D             = tempdata;
        S.fsample_new   = GLEAN.timeseries.settings.fsample;
        S.logtrans      = GLEAN.timeseries.settings.log;
        if isfield(GLEAN.timeseries.settings,'freqbands')
            S.freqbands = GLEAN.timeseries.settings.freqbands;
        else
            S.freqbands = [];
        end
        S.demean    = 0;
        S.prefix    = 'h';

        if strcmp(GLEAN.timeseries.settings.method, 'hilbenv')
            S.method = 'hilbenv';
        elseif strcmp(GLEAN.timeseries.settings.method,'raw')
            S.method = 'none';
        else
            error('GLEAN.timeseries.settings.method not recognised.');
        end
        D = glean_ts_compute(S);

        if GLEAN.timeseries.settings.isepoched
            objpath = [D.path '/C' D.fname];
            D = stacktrialsmeeg(D,objpath,GLEAN.timeseries.settings.conditions);
        end

        % Rename file
        move(D,GLEAN.timeseries.data{session});

        % Remove temporary file
        system(['rm ' tempdata '.*at'])

    end

end


