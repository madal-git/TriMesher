function [X,Y,N_H,N_L,N_I,FACE,O,NODE,N,CELL,M,h,X1,Y1,X2,Y2,FM]=IRT_main(N_H,L,RatioLH,S,D_cell,D_node,D_face)
% [X,Y,N_H,N_L,N_I,NODE,N,CELL,M,h,X1,Y1,X2,Y2,FM]=IRT_main(N_H,L,RatioLH,S,D_cell,D_node,D_face)
% generates a 2D Isosceles Right Triangle mesh on a rectangular domain.
% N_H is the wished number of nodes on height side.
% L is the length of the rectangular domain in the mesh space, not the
% expected value.
% RarioLH is the ratio of length to height.
% S defines the ratio of real dimension to the mesh space dimension
% (e.g. The real length is L*Scale). 
% D_cell is the dimension in the cell data structure.
% D_node is the dimension in the node data structure.
% D_face is the dimension in the face data structure.
% X is x coordinate of all nodes in real space.
% Y is y coordinate of all nodes in real space.
% N_H is the total number of nodes on vertical edge.
% N_L is the total number of nodes on horizontal edge.
% N_I is the total number of interior nodes.
% FACE is the face data structure.
% O is the total number of faces.
% NODE is the node data structure.
% N is the total number of nodes.
% CELL is the cells in real space.
% M is the total number of cells.
% h is the node spacing.
% X1 is the position of left surface.
% X2 is the position of right surface.
% Y1 is the position of bottom surface.
% Y2 is the position of the top surface.
% FM is the flag for mesh type. FM=0 for IRT mesh, this flag is imprtant
% for FVDBM solver.

% Flag for IRT mesh
FM=0;
% Check the ratio input
if RatioLH<1 || RatioLH>80
    error('RatioLH must be within the range between 1 and 80!');
end

if rem(RatioLH,2)~=0 && rem(RatioLH,2)~=1
    %error('RatioLH must be an integer');
end
% Check the length input
if L<2 || L>8 
    error('L must be within the range between 2 and 8!');
end

if L>(6/59)*RatioLH+(112/59) % This formula is based on hand calculation.
    error('Your L is too big according to your RatioLH, try a smaller L!');
end

[C,X1,X2,Y1,Y2]=ini_tri(L,RatioLH,D_cell);  % Initial triangle is added. The cell data structure is defined and the rectangular domain is fixed.

[X_B,Y_B,N_L,N_B,h]=B_node(N_H,L,RatioLH,X1,X2,Y1,Y2); % Generate nodes on the boundaries
Message=['Total ', num2str(N_B), ' nodes are created on the boundaries.'];
disp(Message);

[X_I_B,Y_I_B,N_I_B]=I_nodeB(N_H,N_L,X1,Y1,h); % Generate the interior nodes that are aligned with boundary nodes.
[X_I_AB,Y_I_AB,N_I_AB]=I_nodeAB(N_H,N_L,X1,Y1,h,N_I_B); % Generates the interior nodes that are not aligned with boundary nodes.
N_I=N_I_B+N_I_AB;
Message=['Total ', num2str(N_I), ' nodes are created internally.'];
disp(Message);

[X,Y,N]=node_combine(X_B,Y_B,N_B,X_I_B,Y_I_B,N_I_B,X_I_AB,Y_I_AB,N_I_AB); % Put the coordinates of all nodes into X and Y.
Message=['Total ', num2str(N), ' nodes are created.'];
disp(Message);
%% Generate the cell data structure
[C,M]=IRTtriangulation(N_H,N_L,C,X,Y,N,X_I_AB,Y_I_AB,h,X1,Y1,X2,Y2,D_cell); % Triangulation is done.
[CELL,M]=extract(C,M,D_cell,X1,Y1,X2,Y2); % Get rid of initial triangle from the final triangulation.
[X,Y,X_B,Y_B,X_I_B,Y_I_B,X_I_AB,Y_I_AB,CELL,h,X1,Y1,X2,Y2]=scale(X,Y,X_B,Y_B,X_I_B,Y_I_B,X_I_AB,Y_I_AB,CELL,M,h,X1,Y1,X2,Y2,S); % Scaling
Message=['Total ', num2str(M), ' cells are created.'];
disp(Message);
display('Triangulation is finished!')
display('Cell data structure is generated!')
%% Generate the node data structure
NODE=nddata(CELL,M,X,Y,N_L,N_H,N,N_I,0,h,0,0,FM,D_node);
display('Node data structure is generated!')
%% Generate the face date structure
[FACE,O]=fcdata(CELL,NODE,X,Y,X_B,Y_B,N_B,X_I_B,Y_I_B,N_I_B,X_I_AB,Y_I_AB,N_I_AB,N_L,N_H,X1,X2,Y1,Y2,h,FM,D_face);
Message=['Total ', num2str(O), ' faces are created.'];
disp(Message);
display('Face data structure is generated!')
%% Fill the correlated data in CELL, NODE and FACE
[CELL,NODE,FACE]=datafill(N_I,N_L,N_H,X,Y,CELL,NODE,FACE,FM);
display('Correlated data is filled!')
%% Generate the neighbor information
[CELL,FACE]=IRTngb(CELL,NODE,FACE);
display('Neighbor information is generated!')
%% Checking data
display('Checking mesh...')
meshcheck(CELL,NODE,FACE,N_H,N_L,X,Y,FM);
%% save data
save('IRTmesh.mat','FM','CELL','M','NODE','N','FACE','O','N_I','N_H','N_L','h','X','Y','X1','X2','Y1','Y2');
msg='Please rename the file of generated data! Otherwise it will be overwritten!';
disp(msg);

display('IRT meshing is completed, plotting result...')
figure(3);
for r=1:M;
    P=CELL{r};
    plot(P{22},P{23},P{24},P{25},P{26},P{27});
    hold on
end;
axis equal tight
%%%%% Plot cells
for s=1:N;
    plot(X(s),Y(s));
    hold on
end;
