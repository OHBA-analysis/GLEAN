function glean_convert2spm(filename,data,fs)
% Convert [channels x samples] data at sampling rate fs to SPM12
% format, saving the new .mat and .dat files as "filename.*at"

[pathstr,filestr] = fileparts(filename);
if isempty(pathstr)
    pathstr = pwd;
end
filename = fullfile(pathstr,filestr);

blank = fullfile(fileparts(mfilename('fullpath')),'data','blank_meeg');
D = spm_eeg_load(blank);
D = clone(D,filename,[size(data,1),size(data,2),size(data,3)]);
D(:,:,:) = data;
D = fsample(D,fs);
D.save;

end