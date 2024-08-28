function [L,U_n]=en(P1,P2)
% [L,U_n]=en(P1,P2) calculate the length and the unit normal vector of a line
% connecting two points P1 and P2.
% L is the length of the line.
% U_n is the resulting unit normal, row vector
% P1 is the coordinates of point 1, column vector
% P2 is the coordinates of point 2, column vector

% Hint for determing the direction of normal of a edge with two end points 1
% and 2. The right hand rule determines that the edge normal goes to the
% left-hand side direction of the vector from point 1 to point 2.

L=dis(P1,P2);
U_n(1,1)=-(P2(2,1)-P1(2,1))/L;
U_n(1,2)=(P2(1,1)-P1(1,1))/L;