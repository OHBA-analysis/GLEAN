function glean_envelope(GLEAN)
% Runs the envelope stage of GLEAN.
%
% GLEAN_ENVELOPE(GLEAN)
%
% Adam Baker 2015


pretty_string('RUNNING ENVELOPE STAGE')

% Create temporary directory
[~,tmpdir] = fileparts(tempname);
tmpdir = fullfile(GLEAN.envelope.dir,tmpdir);
mkdir(tmpdir);
c = onCleanup(@() system(['rm -r ' tmpdir]));

for session = 1:numel(GLEAN.data)

    % Check if envelope file exists and whether or not to overwrite
    file_exists = exist(GLEAN.envelope.data{session},'file') == 2;
    overwrite   = GLEAN.envelope.settings.overwrite == 1;

    subspace_method = char(intersect(fieldnames(GLEAN.subspace.settings),{'pca','parcellation','voxel'}));
    if strcmpi(subspace_method,'parcellation')
            msg = 'Subspace method is parcellation: skipping voxelwise envelope computation';
            run_stage = false;
    elseif file_exists
        if overwrite
            msg = ['Overwriting existing envelope file: \n' GLEAN.envelope.data{session} '\n'];
            run_stage = true;
        else
            msg = ['Using existing envelope file: \n' GLEAN.envelope.data{session} '\n'];
            run_stage = false;
        end
    else
        msg = ['Creating new envelope file: \n' GLEAN.envelope.data{session} '\n'];
        run_stage = true;
    end
    fprintf(msg);


    if run_stage

        % Make a temporary filename to copy raw data to
        [~,tempdata] = fileparts(tempname);
        tempdata = fullfile(tmpdir,tempdata);

        % Copy data to temporary filename
        copymeeg(GLEAN.data{session},tempdata)

        % Compute envelopes
        S               = [];
        S.D             = tempdata;
        S.fsample_new   = GLEAN.envelope.settings.fsample;
        S.logtrans      = GLEAN.envelope.settings.log;
        if isfield(GLEAN.envelope.settings,'freqbands')
            S.freqbands = GLEAN.envelope.settings.freqbands;
        else
            S.freqbands = [];
        end
        S.demean    = 0;
        S.prefix    = 'h';
        D = glean_hilbenv(S);

        if GLEAN.envelope.settings.isepoched
            objpath = [D.path '/C' D.fname];
            D = stacktrialsmeeg(D,objpath,GLEAN.envelope.settings.conditions);
        end

        % Rename file
        move(D,GLEAN.envelope.data{session});

        % Remove temporary file
        system(['rm ' tempdata '.*at'])

    end

end


