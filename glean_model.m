function glean_model(GLEAN)
% Runs the model stage of GLEAN.
%
% GLEAN = GLEAN_MODEL(GLEAN)
%
% Adam Baker 2015


pretty_string('RUNNING MODEL STAGE')

% Check if envelope file exists and whether or not to overwrite
file_exists = exist(GLEAN.model.model,'file') == 2;
overwrite   = GLEAN.model.settings.overwrite == 1;
if file_exists
    if overwrite
        msg = ['Overwriting model in: \n' GLEAN.model.model '\n'];
        run_stage = true;
    else
        msg = ['Using existing model in: \n' GLEAN.model.model '\n'];
        run_stage = false;
    end
else
    msg = ['Saving model in: \n' GLEAN.model.model '\n'];
    run_stage = true;
end
fprintf(msg);


if run_stage

    % Concatenate data:
    [dataConcat,subIndx,trlIndex,cndIndex,sessle] = glean_concatenate(GLEAN,'concatenate'); %#ok

    % extract ROIs if requested
    if ~isempty(GLEAN.model.settings.hmm.ROI);
        dataConcat=dataConcat(GLEAN.model.settings.hmm.ROI,:);
    end

    % rb: this does the sign-matching across different sessions
    % TO DO: do the cross-subject task-based sign-matching to allow
    % averaing task data
    if GLEAN.model.settings.hmm.sign_matching
        display(['Running sign-disambiguation session wise']);
        sessle=cell2mat(sessle);
        [flips,scorepath,covmats_unflipped] = findflip(dataConcat',sessle,[]);
    end
    
    dataConcat_fl = flipdata(dataConcat',sessle,flips);
    dataConcat=dataConcat_fl';
    clear dataConcat_fl;
    
    
    dataConcat = normalise(dataConcat,2); % TODO: maybe add an option for this

    switch char(intersect(lower(fieldnames(GLEAN.model.settings)),{'hmm','ica'}));

        case 'hmm'
            options = GLEAN.model.settings.hmm.hmmOptions;

            if ~isempty(trlIndex) && GLEAN.model.settings.hmm.trialwise == 1;
                % tell the hmm that we want to work trialwise
                %T = ones(1,size(cndIndex,2)) * sum(trlIndex == 1);
                % here we estimate trials length per trial
                for trials=1:length(cndIndex)
                    T(trials)=sum(trlIndex == trials)
                end
            else
                % One big epoch please
                T = length(trlIndex);
            end

            % Run model
            hmm = glean_infer_hmm(dataConcat,options,T); %#ok
            save(GLEAN.model.model,'hmm','subIndx')
        case 'ica'
            nICs = GLEAN.model.settings.ica.order;
            [ica.tICs,ica.SM,~] = fastica(dataConcat,     ...
                                          'g','tanh',     ...
                                          'lastEig',nICs, ...
                                          'numOfIC',nICs, ...
                                          'approach','symm');
            save(GLEAN.model.model,'ica','subIndx')
    end

end

end
