
function [A] = SlimProcessSTFiles

%store the old directory
oldDir = pwd;

%format specifier for the files to be read in 
FORMATSPECLOAD = ['%*f%*f%*f%f%*f%*f%*f' repmat('%*s%*f%*f%*f%*f%*f',[1,3]) ...
    '%f%f%s%f\r\n'];

HDR = 'OriginFileName\tImage Name\tX-TP Transformed (pixels, 0,0 is image centre)\tY-TP Transformed (pixels, 0,0 is image centre)\r\n';
FORMATSPECSAVE = '%s\t%s\t%f\t%f\r\n';
SAVEFNM = 'Overall_Results.txt';

%select the files
[fNm fDir] = uigetfile('*Results_*.txt','MultiSelect','On');

%change to the file directory
cd(fDir);

%determine the number of demographic files and ensure file names are stored in a cell array
if iscell(fNm)
    numResultsFiles = length(fNm);
else
    numResultsFiles = 1;
    fNm = {fNm};
end

%output variable 'A'.
A=[];

%load up the selected results files
for i = 1:numResultsFiles
    
    %load in the next file
    fid = fopen(fNm{i},'r');
    tempData = textscan(fid,FORMATSPECLOAD,'Delimiter','\t','HeaderLines',1);
    fclose(fid);
    
    %determine the total number of trials in the data set.
    trialTotInFile = length(tempData{1});
    
    %1: results file name including demographic file name; 2: transformed
    %x-TP; 3: transformed y-TP; 4: image file name (without the extension)
    tempCell = cell(4,trialTotInFile);
    tempFNm = fNm{i}(7:end-4);
    
    %the number of valid trials
    validCount=0;
    
    %populate 'tempCell' (used to append to the cell that is ultimately saved to disk)
    for j=1:trialTotInFile
        
        %store non-bad-flag trials in which the touchpoint does not fall outside the mask
        if tempData{1}(j)==0 && tempData{5}(j)==0
            validCount = validCount+1;
            tempCell{1,validCount} = tempFNm;
            tempCell{2,validCount} = tempData{4}{j,:};
            tempCell{3,validCount} = tempData{2}(j);
            tempCell{4,validCount} = tempData{3}(j);         
        end
    end
    
    %trim 'tempCell'
    if validCount < j
        tempCell(:,validCount+1:end) = [];
    end
    
    %append 'tempCell' to 'A'
    A = [A tempCell];
end

%save 'A' to disk

fid = fopen(SAVEFNM,'w');
fprintf(fid, HDR);
fprintf(fid, FORMATSPECSAVE, A{:,:});
fclose(fid);

%return to the original directory
cd(oldDir);
return;