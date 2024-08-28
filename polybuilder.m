f = msgbox('Please select the Mesher folder location');
uiwait(f);
desiredFolderloc = uigetdir;
desiredFolder = strcat(desiredFolderloc);
currentFolder = strcat(pwd);
if strcmp(currentFolder,desiredFolder) == 0
 cd(desiredFolderloc);
end

dN = Nd';
geS = Seg';

% if Num_holes == 2 || Num_holes == 3
%     eloH = Hole';
% end
if Num_holes ~= 0 
    eloH = Hole';
end
filename = 'PreTri_NewPoly.poly';

fid = fopen(filename,'w');
[L,c] = size(Nd);
addon1 = ' 2 0 1 ';
addon2 = ' 1 ';

formatSpec1 = '%g %g %g %g\n';
formatSpec2 = '%g %g %g \n';

fprintf(fid, num2str(L));
fprintf(fid,addon1);
fprintf(fid, '\n');
fprintf(fid, formatSpec1, dN);

fprintf(fid, num2str(L));
fprintf(fid,addon2);
fprintf(fid, '\n');
fprintf(fid, formatSpec1, geS);

if Num_holes ~= 0
    fprintf(fid, num2str(Num_holes));
    fprintf(fid, '\n');
    fprintf(fid, formatSpec2, eloH);
else
    fprintf(fid, num2str(Num_holes));
    fprintf(fid, '\n');
    fprintf(fid, formatSpec2, Hole);
end

copyfile 'PreTri_NewPoly.poly' Triangle/;

%{

inputsample = false;

while inputsample == false
    prompt = ['Please enter maximum triangle area constraint for mesher' newline '    Value should be between 0 and 1 (default is 0.1)']; 
    maxarea = inputdlg(prompt,'Max Area Constraint');
    maxarea = char(maxarea);
    if(str2num(maxarea) > 0 && str2num(maxarea)< 1)
        inputsample = true;
    else 
        inputsample = false;
    end
end
%}        

seglength = H/(N_H-1);
maxarea = (4.5*seglength^2/5);
maxarea = num2str(maxarea)
cd Triangle/;
command = strcat('triangle -q33a',maxarea,'DenYV PreTri_NewPoly.poly');
status = system(command);
pause(1);

cd ..;
