function GLEAN = glean_settings(GLEAN)
% Check GLEAN options and set defaults.
%
% GLEAN = GLEAN_SETTINGS(GLEAN)
%
% -------------------------------------------------------------------------
%                           ENVELOPE SETTINGS
% -------------------------------------------------------------------------
% .overwrite    - Force overwrite existing files  
%                   permitted: 0/1
%                     default: 0 
% .log          - Apply log transform to envelope 
%                   permitted: 0/1
%                     default: 0 
% .fsample      - New sampling rate of envelope in Hz 
%                   permitted: real number
%                     default: 10
% .freq_bands   - Multiple frequency bands to compute envelope within 
%                   permitted: cell array of [Hz_low Hz_high]
%                     default: {[0 Inf]}
%
% -------------------------------------------------------------------------
%                           SUBSPACE SETTINGS
% -------------------------------------------------------------------------
% .overwrite        - Force overwrite existing files 
%                       permitted: 0/1
%                         default: 0 
% .normalisation    - Type of normalisation to apply (before computing the 
%                     subspace. 'none': no normalisation applied. 'voxel': 
%                     normalise each voxel to have a standard deviation of 
%                     one. 'global': normalise each session to have an 
%                     average standard deviation of one over all voxels.
%                       permitted: 'none','voxel','global'
%                         default: 'none'
% 
% Plus one of the following subspace methods with additional options:
% .pca
%   .dimensionality - Number of principal components to retain
%                       permitted: Integer less than the rank of the data 
%                         default: 40 
%   .whiten         - Remove variance of principal components
%                       permitted: 0/1
%                         default: 0 
% .parcellation
%   .file               - A file .nii or .nii.gz file at the same gridstep 
%                         as the data OR a .mat file containing a matrix of 
%                         dimensions voxels x parcels with the same number 
%                         of voxels as the data.
%                           permitted: string containing a valid file
%                           default: No default (must be specified)
%   .mask               - If the parcellation is in nifti format, then also
%                         need to specift a .nii or .nii.gz file containing 
%                         a wholebrain mask.
%                           permitted: string containing a valid file
%                           default: No default
%   .orthogonalisation  - Type of orthogonalisation to apply (prior to 
%                         computing parcelwise envelopes). See help text in 
%                         remove_source_leakage.m
%                           permitted: 'none','symmetric','closest',
%                                      'householder'
%                             default: 'none'
%   .method             - Type of method to use for reconstructing parcel
%                         time courses. See help text in get_node_tcs.m
%                           permitted: 'pca','mean','peakvoxel',
%                                      'spatialbasis'
%                             default: 'peakvoxel'
%
% -------------------------------------------------------------------------
%                             MODEL SETTINGS
% -------------------------------------------------------------------------
% .overwrite - Force overwrite existing files  
%                permitted: 0/1
%                  default: 0
% Plus one of the following model types with additional options:
% .hmm
%   .nstates - Number of states to infer
%                permitted: Integer greater than 1
%                  default: 8
%   .nreps   - Number of times to run the inference (model with the 
%              greatest model evidence is retained.
%                permitted: Integer greater than 0
%                  default: 1
% .ica       
%   .order   - Number of components to retain
%                permitted: Integer greater than 0
%                  default: 20
%
% -------------------------------------------------------------------------
%                            RESULTS SETTINGS
% -------------------------------------------------------------------------
% One of the following result types with additional options:
% .pcorr OR .connectivity_profile
%    .format       - Format to save maps as
%                      permitted: 'mat','nii'
%                        default: 'mat'
%    .space        - Subspace to save maps as
%                      permitted: 'parcel','voxel',{'parcel','voxel'}
%                        default: 'voxel'
%
% -------------------------------------------------------------------------
%
% Adam Baker 2015

% Loop through each module and set up the directory and file structure
for module = {'envelope','subspace','model','results'}
    % Check if module exists and initialise if not
    if ~isfield(GLEAN,char(module)) || ~isstruct(GLEAN.(char(module)))
        GLEAN.(char(module)) = struct;
    end
    GLEAN = check_settings(GLEAN,module);
end

end


function GLEAN = check_settings(GLEAN,module)

if ~isfield(GLEAN.(char(module)),'settings') || ~isstruct(GLEAN.(char(module)).settings)
    GLEAN.(char(module)).settings = struct;
end

switch char(module)
    
    case 'envelope'
        % --- VALIDATE ENVELOPE FIELDS --- %            
        V = {};
        V = addOption(V,'overwrite', ...
                        0, ...
                        0, ...
                        @(x) x==0 || x == 1);
        
        V = addOption(V,'log', ...
                        0, ...
                        0, ...
                        @(x) x==0 || x == 1);

        V = addOption(V,'fsample', ...
                        0, ...
                        10, ...
                        @(x) isnumeric(x) && isscalar(x) && (x > 0));

        V = addOption(V,'freqbands', ...
                        0, ...
                        {[0 Inf]}, ...
                        @(x) iscell(x) && all(cellfun(@isnumeric,x)));
                   
        GLEAN.envelope.settings = validateOptions(V,GLEAN.envelope.settings,'GLEAN.envelope.settings'); 
   
    case 'subspace'
        % --- VALIDATE COMMON SUBSPACE FIELDS --- %
        subspace = intersect(fieldnames(GLEAN.subspace.settings),{'pca','voxel','parcellation'});
        if numel(subspace) ~= 1 || ~isstruct(GLEAN.subspace.settings.(char(subspace)))
            error('GLEAN.subspace.settings should contain one and only one of the following structures: ''pca'', ''voxel'', or ''parcellation''');
        else
            subspace = char(subspace);
        end      
                     
        V = {};
        V = addOption(V,'overwrite', ...
                      0, ...
                      0, ...
                      @(x) x==0 || x == 1);
                   
        V = addOption(V,'normalisation', ...
                      0, ...
                      'none', ...
                      @(x) any(strcmpi(x,{'none','voxel','global'})));

        V = addOption(V,subspace, ...
                        0, ...
                        struct, ...
                        @isstruct);
                              
        GLEAN.subspace.settings = validateOptions(V,GLEAN.subspace.settings,'GLEAN.subspace.settings'); 


        % --- VALIDATE SPECIFIC SUBSPACE FIELDS --- %
        switch subspace
            case 'pca'
                V = {};
                V = addOption(V,'dimensionality', ...
                              0, ...
                              40, ...
                              @(x) isnumeric(x) && isscalar(x) && (x > 0));

                V = addOption(V,'whiten', ...
                              0, ...
                              0, ...
                              @(x) x==0 || x == 1);
                GLEAN.subspace.settings.(subspace) = validateOptions(V,GLEAN.subspace.settings.(subspace),['GLEAN.subspace.settings.' subspace]);

            case 'parcellation'
                V = {};
                V = addOption(V,'file', ...
                              1, ...
                              [], ...
                              @(x) ischar(x) || iscellstr(x));

                V = addOption(V,'mask', ... 
                              0, ...
                              [], ...
                              @(x) ischar(x) || iscellstr(x));
                          
                V = addOption(V,'method', ...
                              0, ...
                              'peakvoxel', ...
                              @(x) any(strcmpi(x,{'pca','mean','peakvoxel','spatialbasis'})));
           
                V = addOption(V,'orthogonalisation', ...
                              0, ...
                              'none', ...
                              @(x) any(strcmpi(x,{'none','symmetric','closest','householder'})));
                
                    GLEAN.subspace.settings.(subspace) = validateOptions(V,GLEAN.subspace.settings.(subspace),['GLEAN.subspace.settings.' subspace]);


            case 'voxel'
                GLEAN.subspace.settings.(subspace) = [];

        end


        
    case 'model'
        % --- VALIDATE COMMON MODEL FIELDS --- %
        model = intersect(fieldnames(GLEAN.model.settings),{'hmm','ica'});
        if numel(model) ~= 1 || ~isstruct(GLEAN.model.settings.(char(model)))
            error('GLEAN.model.settings should contain one and only one of the following structures: ''ica'', ''hmm''');
        else
            model = char(model);
        end    
        V = {};
        V = addOption(V,'overwrite', ...
                      0, ...
                      0, ...
                      @(x) x==0 || x == 1);
                   
        V = addOption(V,model, ...
                        0, ...
                        struct, ...
                        @isstruct);
                              
        GLEAN.model.settings = validateOptions(V,GLEAN.model.settings,'GLEAN.model.settings'); 
        
        % --- VALIDATE SPECIFIC MODEL FIELDS --- %                
        V = {};
        switch model
            case 'hmm'
                V = addOption(V,'nstates', ...
                                0, ...
                                8, ...
                                @(x) isnumeric(x) && isscalar(x) && (x > 0));
                
                V = addOption(V,'nreps', ...
                                0, ...
                                1, ...
                                @(x) isnumeric(x) && isscalar(x) && (x > 0));
                
            case 'ica'
                V = addOption(V,'order', ...
                                0, ...
                                20, ...
                                @(x) isnumeric(x) && isscalar(x) && (x > 0));
                
        end
        GLEAN.model.settings.(model) = validateOptions(V,GLEAN.model.settings.(model),['GLEAN.model.settings.' model]);



    case 'results'
        % --- VALIDATE COMMON RESULTS FIELDS --- %
        result_types = intersect(fieldnames(GLEAN.results.settings),{'pcorr','connectivity_profile'});

        for result_type = result_types
            result = char(result_type);
            switch result
                case {'pcorr','connectivity_profile'}
                    V = {};
                    V = addOption(V,'format', ...
                                    0, ...
                                    'mat', ...
                                    @(x) any(strcmpi(x,{'mat','nii'})));
                          
                    V = addOption(V,'space', ...
                                    0, ...
                                    'voxel', ...
                                    @(x) any(strcmpi(x,{'voxel','parcel'})));
                    GLEAN.results.settings.(result) = validateOptions(V,GLEAN.results.settings.(result),['GLEAN.results.settings.' result]);
      
            end
                               
        end
              
end



end



function S = addOption(S,option,required,default,validFcn)

i = numel(S) + 1;
S(i).option     = option;
S(i).required   = required;
S(i).default    = default;
S(i).validFcn   = validFcn;

end

 
function S = validateOptions(V,S,callerID)
% Checks input options. Throws an error if any extra options are specified 
% and notify if any default options are used. Matlab has inbuilt 
% functionality for this but it is very buggy.
    
    vOpts = {V.option};
    sOpts = fieldnames(S);
    
    % Check each option in S for invalid options:
    for i = 1:length(sOpts)
        v = V(strcmp({V.option},sOpts(i))); % select validator
        if isempty(v)
            % Check if option is present in V:
            error(['Unrecognised option ''' char(sOpts(i)) ''' in ' callerID]);
        elseif v.validFcn(S.(char(sOpts(i)))) == 0
            % Check if option is present in V:
            error(['Option ''' v.option ''' in ' callerID ' must satisfy: ' func2str(v.validFcn)]);
        end
    end
        
    % Check each option in V for undefined options:  
    defaultStr = '';
    for i = 1:length(vOpts)
        v = V(i); % select validator
        if ~ismember(v.option,sOpts)
            if v.required % Throw error if input is required
                error(['Must specify ' v.option ' in ' callerID])
            else % Set default if option is not defined
                S.(v.option) = v.default;
                if ~isempty(defaultStr)
                    defaultStr = [defaultStr ', ']; %#ok
                end
                defaultStr = [defaultStr v.option]; %#ok
            end
        end
    end
    if ~isempty(defaultStr)
        disp(['Using default values for ' defaultStr ' in ' callerID]);
    end
    
end





