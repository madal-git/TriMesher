%If code is run as part of PreTri_poly_main, will use GT
%If code is run as a standalone program, will use IRT

tester = exist('GT','var');

if(tester ~= 0)
    fm = 1;
else
    fm = 0;
end

%fm=1; % Flag for mesh type. fm=0---IRT mesh; fm=1---Triangle mesh
IRT_N_H=36;
IRT_L=2;
IRT_RatioLH=1;
IRT_S=1;
% IRT_N_H is the wished number of nodes on height side.
% IRT_L is the length of the rectangular domain in the mesh space, not the
% expected value.
% IRT_RarioLH is the ratio of length to height.
% IRT_S defines the ratio of real dimension to the mesh space dimension
tic;
if fm==0
    [X,Y,N_H,N_L,N_I,FACE,O,NODE,N,CELL,M,h,X1,Y1,X2,Y2,FM]=IRT_main(IRT_N_H,IRT_L,IRT_RatioLH,IRT_S,40,25,55);
elseif fm==1
    fileloads;
    msg1=['1. Copy the mesh file .ele, .node, .neigh and .edge into the current working directory.'];
    disp(msg1);
    msg2=['2. Create four .txt files in the current working directory named as ele, node, neigh and edge seperately.'];
    disp(msg2);
    msg3=['3. Open .ele file with text editor and copy its content into ele.txt file, then delete the first and last line then save.'];
    disp(msg3);
    msg4=['4. Repeat for node, neigh and edge files.'];
    disp(msg4);
    msg5=['5. Import data from four .txt files into Matlab. Select Delimited and and select Numeric Matrix for Output Type.'];
    disp(msg5);
    msg6=['6. Save the imported raw mesh data as mat file for future use.'];
    disp(msg6);
    save('ImportedTri.mat','ele','node','neigh','edge');
    msg='Please rename the file of PostTri data! Otherwise it will be overwritten when new data is imported again!';
    [FM,CELL,M,NODE,N,FACE,O,N_I,N_I_N,N_H,N_L,dx,dy,X,Y,X1,X2,Y1,Y2]=PostTri_rec_main(ele,node,neigh,edge);
    %% delete raw mesh data in the memory and other temperary variables
    clear('msg1');
    clear('msg2');
    clear('msg3');
    clear('msg4');
    clear('msg5');
    clear('msg6');
    clear('msg');
    clear('ele');
    clear('node');
    clear('neigh');
    clear('edge');
else
    error('Wrong flag for mesh type!');
end
toc;
%% delete raw mesh data in the memory and other temperary variables
clear('fm');
clear('IRT_N_H');
clear('IRT_L');
clear('IRT_RatioLH');
clear('IRT_S');