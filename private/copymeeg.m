function copymeeg(infile,outfile)

% Copy data to temporary filename
[inpath,infile] = fileparts(infile);
[outpath,outfile] = fileparts(outfile);

for ext = {'.mat','.dat'}
    system(['cp ' fullfile(inpath,infile) char(ext) ' ' fullfile(outpath,outfile) char(ext)]);
end

load(fullfile(outpath,outfile)); % loads D

D.path = outpath;
D.fname = [outfile,'.mat'];
D.data.fname = [fullfile(outpath,outfile),'.dat'];

save(fullfile(outpath,outfile),'D');

D = spm_eeg_load(fullfile(outpath,outfile));
if ~islinked(D)
    error('MEEG object is not linked to a valid datafile')
end

end