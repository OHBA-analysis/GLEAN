% Load GLEAN
load('/Users/abaker/Scratch/GLEAN/Schizophrenia/glean.mat');

TBL = readtable('/Users/abaker/Data/Schizophrenia/demographics/data-proc-notes.xlsx');

ControlID = TBL(TBL.Include==1 & strcmp(TBL.Disorder,'Control'),:).Participant;
PatientID = TBL(TBL.Include==1 & strcmp(TBL.Disorder,'Schizophrenia'),:).Participant;


sessionIDs = regexp(GLEAN.data,'\d\d\d\d','match');
sessionIDs = str2double([sessionIDs{:}]'); 

group_labels = cell(size(GLEAN.data));
group_labels(ismember(sessionIDs,ControlID)) = {'Control'};
group_labels(ismember(sessionIDs,PatientID)) = {'Patient'};

groups = unique(group_labels)';

% Unpaired test design matrix
design = zeros(length(group_labels),length(groups));
for group = groups
    design(ismember(group_labels,group),ismember(groups,group)) = 1;
end

contrast = [-1 1]; % Patients - Controls






%% Group differences in mean of envelopes:
M = zeros(num_channels,num_subs);
for sub = 1:num_subs
    
    % Get weights normalisation:
    D = spm_eeg_load(GLEAN.data{sub});
    montage1 = montage(D,'getmontage',1);
    montage2 = montage(D,'getmontage',2);
    weights_norm = montage2.tra(:,1)./montage1.tra(:,1);
    
    % Mean of envelope (M*f) 
    D = spm_eeg_load(GLEAN.envelope.data{sub});
    M(:,sub) = mean(D(:,:,:,:),3); % TODO: maybe block this for memory
    
    % Remove weights normalisation:
    M(:,sub) = M(:,sub) ./ repmat(weights_norm,1,D.nfrequencies);
    
end

input_nii = fullfile(savedir,'mean_envelope.nii.gz');
nii_quicksave(M,input_nii,8)

output_nii = fullfile(savedir,'randomise');


command = sprintf('randomise -i %s -o %s -d %s -t %s -T', ...
                  input_nii, output_nii, design_file, contrast_file);

system(command)
              
%randomise -i nii_file -o TwoSampT -d design.mat -t design.con -m mask -T 
%randomise -i <4D_input_data> -o <output_rootname> -d design.mat -t design.con -m <mask_image> -n 500 -D -T
%tmp=['randomise -d ' dirname 'design.mat -t ' dirname 'design.con -i ' fname_cope ' -o ' output_name ' -1 -T -R -n ' num2str(nperms) ' --seed=0 -v ' num2str(group_varcope_spatial_smooth_fwhm) ' -m ' mask_fname]; % -c means cluster-based thresholding


%copes - [vox x subs]


% Fit GLM:
[~,~,~,tstat] = glean_glm(M',design,contrast);
nii_quickview(tstat,8)




%% Group differences in source variance:
V = zeros(num_subs,num_channels);
for sub = 1:num_subs
    D = spm_eeg_load(GLEAN.data{sub});
    D = montage(D,'switch',1);
    V(sub,:) = osl_source_variance(D);
end

% Fit GLM:
[~,~,~,tstat] = glean_glm(V,design,contrast);
nii_quickview(tstat',8)





%% Group level differences in HMM state COPEs
model = load(GLEAN.model.model);
copes = zeros(num_subs,num_channels,model.hmm.K);
for sub = 1:num_subs
        
    % Get subject statepath:
    sp_sub = cell2mat(arrayfun(@(k) ...
                      model.hmm.statepath(model.subIndx == sub) == k, ...
                      1:model.hmm.K, ...
                      'UniformOutput', 0));
    
    % Get weights normalisation:
    D = spm_eeg_load(GLEAN.data{sub});
    montage1 = montage(D,'getmontage',1);
    montage2 = montage(D,'getmontage',2);
    weights_norm = montage2.tra(:,1)./montage1.tra(:,1);
    
    % Load envelope data
    D = spm_eeg_load(GLEAN.envelope.data{sub});
    
    % Weights un-normalised COPEs
    copes(sub,:,:) = glean_regress(D,sp_sub,'cope',1./weights_norm);
    
end

FWE_corrected_tstats = zeros(num_channels,model.hmm.K);
for k = 1:model.hmm.K
    % Save COPE for each state
    input_nii  = fullfile(savedir,sprintf('cope%i.nii.gz',k));
    output_nii = fullfile(savedir,sprintf('randomise%i',k));
    writenii(copes(:,:,k)',input_nii,GLEAN.envelope.settings.mask);
    command = sprintf('randomise -i %s -o %s -d %s -t %s -x', ...
                      input_nii, output_nii, design_file, contrast_file);
    fprintf('Running FSL randomise for state %i of %i\n',k,model.hmm.K)
    [~,~] = system(command);
    FWE_corrected_tstats(:,k) = readnii([output_nii,'_vox_corrp_tstat1.nii'], ...
                                        GLEAN.envelope.settings.mask);
    system(['rm ' output_nii '*.*']);
end
writenii(FWE_corrected_tstats,fullfile(savedir,'FWE_corrected_tstats.nii.gz'),GLEAN.envelope.settings.mask);

              



% Fit GLM for each state:
tstat = zeros(num_channels,model.hmm.K);
for k = 1:model.hmm.K
    [~,~,~,tstat(:,k)] = glean_glm(copes(:,:,k),design,contrast);
end
nii_quickview(tstat,8)


%% Group level differences in HMM state covariance
model = load(GLEAN.model.model);
netmat = zeros(num_subs,num_channels,model.hmm.K);

% Concatenate data:
dataConcat = cell(1,num_subs); % cell arrays grow better than arrays
for sub = 1:num_subs
    D = spm_eeg_load(GLEAN.subspace.data{sub});
    dat = reshape(D(:,:,:,:),[],D.nsamples,D.ntrials);
    dataConcat{sub} = dat;
end
dataConcat = cell2mat(dataConcat);
 
%%
figure
map = zeros(num_channels,size(dataConcat,1));
for k = 1:model.hmm.K;
    netmat = cov(dataConcat(:,model.hmm.statepath==k)');
    map(:,k) = parcellation2map(diag(netmat),GLEAN.subspace.settings.parcellation.file,GLEAN.envelope.settings.mask);
    subplot(2,4,k)
    osl_braingraph(corrcov(netmat), [], diag(netmat), [], mni_coords, [], 99)
end
%nii_quickview(map,8)

%%

parcellation = nii_quickread(GLEAN.subspace.settings.parcellation.file,8);
mni_coords = osl_mnimask2mnicoords(GLEAN.envelope.settings.mask);
[~,i] = max(parcellation);
mni_coords = mni_coords(i,:);

 

















