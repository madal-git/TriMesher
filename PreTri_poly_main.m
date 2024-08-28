H=2;
RatioL2H=2/H; %For 1:1 ratio set top number to be same as H
N_H=30; % Nodes in height
N_L=30; % Nodes in length
Origin=[0;0];
Fd=2; % Flag for adjusting the node spacing on the outer rectangular boundaries
FBuf_inner=3; % Flag for buffer laylers on the inner rectangular boundaries
% This value will be automatically adjusted to smoothen the cell size in the layer and the 
% size in the surrounding area if nonzero 
FBuf_outer=2; % Flag for buffer laylers on the outer rectangular boundaries
Fw=0; % Flag for whether the line of nodes in the wake region is needed. 0---No; 1---Yes
Fcob=0; % Flag for whether the outer boundaries are a customized shape (a non-customized shape is a regular rectangular shape). 0---No (Rectangular); none-zero---Yes
Num_holes=1 ; % Change

hole1=cell(5,1);

hole1{1}=1; % 1--circlular; 0--Rectangular
hole1{2}=H/6; % Radius if circular
hole1{3}=H/2;  % Hole center to the left boundary
hole1{4}=H/2; % Hole center to the bottom boundary
hole1{5}=40; % Number of nodes on hole

hole2=cell(5,1);

hole2{1}=1;
hole2{2}=0.3;
hole2{3}=3;
hole2{4}=4;
hole2{5}=11;
 
% hole3=cell(5,1);
% hole3{1}=1;
% hole3{2}=0.3;
% hole3{3}=1.6;
% hole3{4}=0.4;
% hole3{5}=200;

if Num_holes == 2
holes_para=[hole1,hole2];
elseif Num_holes == 3
holes_para={hole1,hole2,hole3};
else
holes_para=[hole1];
end

[Nd, Seg, Hole]=PreTri_poly_rec(H,RatioL2H,N_H,N_L,Origin,Fd,Num_holes,holes_para,FBuf_inner,FBuf_outer,Fcob,Fw);

%%%% Instructions for preparing the data for Triangle
%{
msg='Please proceed the following steps to preapre the data for Triangle.';
disp(msg);
msg1='1. Create a .txt file with desired name properly describing the nature of the mesh (the name should start with PreTri_.....).';
msg2='2. Open the .txt file.';
msg3='3. Copy and paste the single line for nodes in the result window.';
msg4='4. Press Enter to start a new line. Then open the variable Nd in the Workspace, then copy and paste all entries.';
msg5='5. In the same file, tab a new line and repeat step 3 & 4 for segments (variable named Seg).';
msg6='6. In the same file, tab a new line and repeat step 3 & 4 for holes (variable named Hole).';
msg7='7. Close the text file and change the extension from .txt to .poly and move the file to the working directory of Triangle.';
msg8='8. Open the command window in Linux and navigate to the directory where Triangle is saved.';
msg9='9. Type in the command "triangle -q33a.1DenYV YOUR_FILE_NAME.poly" and pressure enter to mesh. Use the command showme YOUR_FILE_NAME.poly to see the generated mesh.';
msg10='10. For detailed manual of using Triangle, visit https://www.cs.cmu.edu/~quake/triangle.html.';
msg11='11. Change the file names to include "PostTri" for clarfication. Copy the files YOUR_FILE_NAME.1.ele, YOUR_FILE_NAME.1.edge, YOUR_FILE_NAME.1.neigh and YOUR_FILE_NAME.1.node to Matlab working directory for PostTri functions.';

disp(msg1);
disp(msg2);
disp(msg3);
disp(msg4);
disp(msg5);
disp(msg6);
disp(msg7);
disp(msg8);
disp(msg9);
disp(msg10);
disp(msg11);
%}

polybuilder;
GT = 'true';
mesher;