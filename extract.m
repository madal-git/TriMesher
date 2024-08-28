function [CELL,M]=extract(C,M,D,X1,Y1,X2,Y2)
% function [CELL,M]=extract(C,M,D) gets rid of the initial triangle from the
% triangulation.
% C and M are defined in the function IRTtriangulation.
% D is defined in the main function.
% X1,Y1,X2,Y2 are defined in the function ini_tri.
% CELL is the finial triangulation.
% M is the total number of triangles.

CELL=cell(D,1);
CELL={CELL,CELL};
b=0;
for p=1:M;
    T=C{p+1};
    Centroid=T{5};
    if Centroid(1,1)>X1 && Centroid(1,1)<X2 && Centroid(2,1)>Y1 && Centroid(2,1)<Y2;
        b=b+1;
        T{1}=b;
        %%%%% Renumber the cell number
        CELL{b}=T;
    else
        ;
    end;
end;
if b~=M
    error('cell missing!');
end