f = msgbox('Please select the Triangle folder location');
uiwait(f);
desiredFolderloc = uigetdir;
desiredFolder = strcat(desiredFolderloc);
currentFolder = strcat(pwd);
if strcmp(currentFolder,desiredFolder) == 0
 cd(desiredFolderloc);
end

copyfile PreTri_NewPoly.1.edge ../edge.txt;
copyfile PreTri_NewPoly.1.ele ../ele.txt;
copyfile PreTri_NewPoly.1.neigh ../neigh.txt;
copyfile PreTri_NewPoly.1.node ../node.txt;

cd ../;

edgein = importdata('edge.txt', ' ',1);
edge = edgein.data(1:end,:);
elein = importdata('ele.txt', ' ',1);
ele = elein.data(1:end,:);
neighin = importdata('neigh.txt', ' ',1);
neigh = neighin.data(1:end,:);
nodein = importdata('node.txt', ' ',1);
node = nodein.data(1:end,:);

clear vars currentFolder desiredFolder edgein elein neighin nodein;

