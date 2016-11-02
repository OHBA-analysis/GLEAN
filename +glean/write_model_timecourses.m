function write_model_timecourses(GLEAN)
% Write GLEAN model time courses as SPM12 MEEG objects.
%
% glean.write_model_timecourses(GLEAN)
%
% Adam Baker 2015   
    model = load(GLEAN.model.model);
    
    for session = 1:numel(GLEAN.data)
            
        D = spm_eeg_load(GLEAN.subspace.data{session});
        
        if isfield(model,'hmm')
            K = model.hmm.K;
            D = clone(D,GLEAN.model.data{session},[2*K,D.nsamples,D.ntrials]);
            D(1:K,:,:) = get_statepath(model,session,'hard');
            D(K+1:D.nchannels,:,:) = get_statepath(model,session,'soft');
            D = chantype(D,1:K,'Viterbi');
            D = chanlabels(D,1:K,arrayfun(@(k) strcat('state',num2str(k)),1:8,'UniformOutput',0));
            D = chantype(D,K+(1:K),'Gamma');
            D = chanlabels(D,K+(1:K),arrayfun(@(k) strcat('state',num2str(k)),1:8,'UniformOutput',0));
        elseif isfield(model,'ica')
            error('Please implement')
        end
        
        D.save
                    
    end
end


    
function statepath = get_statepath(model,session,type)
   switch(type)
       case 'hard'
           statepath = cell2mat(arrayfun(@(k) ...
                                model.hmm.statepath(model.subIndx == session) == k, ... 
                                1:model.hmm.K, ...
                                'UniformOutput', 0));
       case 'soft'
            statepath = model.hmm.train.Gamma(model.subIndx == session,:);
   end    
   statepath = statepath';
end