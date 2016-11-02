function model(GLEAN)
% Runs the model stage of GLEAN.
%
% GLEAN = glean.model(GLEAN)
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
    [dataConcat,subIndx] = glean.concatenate(GLEAN,'concatenate'); %#ok
    dataConcat = normalise(dataConcat,2); % TODO: maybe add an option for this

    switch char(intersect(lower(fieldnames(GLEAN.model.settings)),{'hmm','ica'}));

        case 'hmm'
            hmmSettings = struct('K',GLEAN.model.settings.hmm.nstates,    ...
                                 'order',0,                               ...
                                 'Ninits',GLEAN.model.settings.hmm.nreps, ...
                                 'zeromean',0);
            hmm = glean.infer_hmm(dataConcat,hmmSettings); %#ok
            save(GLEAN.model.model,'hmm','subIndx')
        case 'ica'
            nICs = GLEAN.model.settings.ica.order;
            [ica.tICs,ica.SM,~] = fastica(dataConcat,            ...
                                          'g','tanh',            ...
                                          'lastEig',nICs,        ...
                                          'numOfIC',nICs,        ...
										  'stabilization', 'on', ...
                                          'approach',      'symm');
            save(GLEAN.model.model,'ica','subIndx')
			
		otherwise
			error([mfilename ':unsupportedModel'], ...
				  'Unrecognised model choice. \n');
    end
    
end

end