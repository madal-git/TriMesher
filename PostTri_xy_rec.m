function [X,Y,X1,X2,Y1,Y2,NH,NL,dx,dy]=PostTri_xy_rec(node,N_I,N)
% function [X,Y,X1,X2,Y1,Y2,NH,NL,dx,dy]=PostTri_xy_rec(node) generates all
% the variables associated with coordinates from node.
% X---the x coordinate of all nodes
% Y---the y coordinate of all nodes
% X1---the position of left surface of outer rectangular domain
% X2---the position of right surface of outer rectangular domain
% Y1---the position of bottom surface of outer rectangular domain
% Y2---the position of top surface of outer rectangular domain
% NH---The number of nodes on left or right surface of outer rectangular domain
% NL---The number of nodes on top or bottom surface of outer rectangular domain
% dx---the node spacing on top or bottom surface of outer rectangular domain
% dy---the node spacing on top or bottom surface of outer rectangular domain

% N_I is the total number of nodes except the outer boundary nodes
% N is the total number of nodes

%% X and Y
for r=1:N
    X(r)=node(r,2);
    Y(r)=node(r,3);
end
%% X1, X2, Y1, Y2
X1=min(X);
X2=max(X);
Y1=min(Y);
Y2=max(Y);
%% NH, NL
% NH
h1=0;
h2=0;
for r=1:N
    if X(r)==X1
        h1=h1+1;
    elseif X(r)==X2
        h2=h2+1;
    else
        ;
    end
end
if h1~=h2
    error('The number of nodes on the left and right surface should be the same!');
else
    NH=h1;
end
% NL
l1=0;
l2=0;
for r=1:N
    if Y(r)==Y1
        l1=l1+1;
    elseif Y(r)==Y2
        l2=l2+1;
    else
        ;
    end
end
if l1~=l2
    error('The number of nodes on the bottom and top surface should be the same!');
else
    NL=l1;
end
% check
if (N-N_I)~=2*NH+2*NL-4
    error('Check the N, N_I, NL and NH!');
end
%% dx, dy
dx=(X2-X1)/(NL-1);
dy=(Y2-Y1)/(NH-1);