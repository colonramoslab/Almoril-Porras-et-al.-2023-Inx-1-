function directoryInfo = getDirectoryInformation(dirname, examineBinFiles)
%function directoryInfo = determineWhichSetup(dirname)

existsAndDefault('examineBinFiles', false);

d = dir(fullfile(dirname, '\*.bin'));
if (isempty(d))
    directoryInfo = [];
    return;
end
directoryInfo.binFileName = '';
directoryInfo.headerName = '';
directoryInfo.timfname = '';
directoryInfo.fstub = '';
directoryInfo.setupNumber = -1;
directoryInfo.error = '';
directoryInfo.x0 = NaN;
directoryInfo.y0 = NaN;
directoryInfo.w = NaN;
directoryInfo.h =NaN;
directoryInfo.metrics = [];

directoryInfo = repmat(directoryInfo, size(d));

for j = 1:length(d)
    directoryInfo(j).binFileName = fullfile(dirname, d(j).name);
    if (strncmp(d(j).name, '2009', 4))
        directoryInfo(j).setupNumber = 1;
    end
    stub = regexp(d(j).name, '.*[^(\.bin)]', 'match');
    directoryInfo(j).headerName = fullfile(dirname,[stub{1} '_header.txt']);
    directoryInfo(j).timfname = fullfile(dirname,[stub{1} '.tim']);
    
end

for j = 1:length(directoryInfo)
    try
        fid = fopen(directoryInfo(j).headerName, 'r');
    catch 
        directoryInfo(j).error = 'could not open file';
        fclose(fid);
        continue;
    end
    if (fid <= 0)
        directoryInfo(j).error = ['fopen returned ' num2str(fid)];
        continue;
    end
    [msg, errnum] = ferror(fid);
    if (errnum ~= 0)
        directoryInfo(j).error = ['fopen had error: ' msg];
        fclose(fid);
        continue;
    end
    
    while (~feof(fid))
        str = fgetl(fid);
        if (~ischar(str))
            break;
        end
        pat = '(fstub:)(.*)';
        tok = regexp(str, pat, 'tokens');
        if (~isempty(tok))
            break;
        end
    end
    if (isempty(tok))
        directoryInfo(j).error = 'could not find "fstub:" in file anywhere';
        fclose(fid);
        continue
    end
    fstub = tok{1}{2};
    directoryInfo(j).fstub = fstub;
   
    issetup1 = ~isempty(regexpi(fstub, '(my computer)', 'once'));
    issetup2 = ~isempty(regexpi(fstub, '(mason)', 'once'));
    
    while(~feof(fid))
        str = fgetl(fid);
        if (~ischar(str))
            break;
        end
        pat = 'analysis rectangle:';
        tok = regexp(str, pat, 'tokens');
        if (~isempty(tok))
            break;
        end
    end
    try 
        directoryInfo(j).x0 = cell2mat(textscan(fgetl(fid), ' - %d'));
        directoryInfo(j).y0 = cell2mat(textscan(fgetl(fid), ' - %d'));
        directoryInfo(j).w = cell2mat(textscan(fgetl(fid), ' - %d'));
        directoryInfo(j).h = cell2mat(textscan(fgetl(fid), ' - %d'));
    catch
        directoryInfo(j).error = [directoryInfo(j).error '  could not find analysis rectangle'];
    end
    
    fclose(fid);
    if (issetup1 && issetup2)
        directoryInfo(j).setupNumber = 3;
        directoryInfo(j).error = '  my computer and mason both in file name';
        continue;
    end
    if (issetup1)
        directoryInfo(j).setupNumber = 1;
        continue;
    end
    if (issetup2)
        directoryInfo(j).setupNumber = 2;
        continue;
    end
    directoryInfo(j).setupNumber = 0;
    directoryInfo(j).error = '  neither mason nor my computer in fstub';
end


for j = 1:length(directoryInfo)
    
    if (directoryInfo(j).setupNumber > 0)
        continue;
    end
    try
        data = importdata(directoryInfo(j).timfname);
    catch 
        directoryInfo(j).error = [directoryInfo(j).error ' could not import data from tim file'];
        continue;
    end
    if (isempty(data))
        directoryInfo(j).error = [directoryInfo(j).error ' trouble reading tim file'];
    end
    dt = median(diff(data(:,1)));
   %{ 
    if (dt > 450 && dt < 550)
        directoryInfo(j).setupNumber = 11;
        directoryInfo(j).error = [directoryInfo(j).error ' from timing, this is a linjiao file'];
    end
    %}
    if (dt > 150 && dt < 350)
        directoryInfo(j).setupNumber = 22;
        directoryInfo(j).error = [directoryInfo(j).error ' from timing, this is a bracher file'];
    end
    
end

if (~examineBinFiles)
    return;
end

matlabpool open

parfor j = 1:length(directoryInfo)
    try
        directoryInfo(j).metrics = getBasicInformationFromBinFile(directoryInfo(j).binFileName);
    catch me
        disp (me.getReport());
        directoryInfo(j).error = [directoryInfo(j).error me.getReport()];
    end
end
    
matlabpool close