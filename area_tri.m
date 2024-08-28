function A=area_tri(nd1,nd2,nd3)
% function A=area_tri(nd1,nd2,nd3) calulate the area of a triangle
% A is the area
% nd1, nd2 and nd3 are the coordinates of three vertices, has to be column
% vector.

A=abs((nd1(1,1)*(nd2(2,1)-nd3(2,1))+nd2(1,1)*(nd3(2,1)-nd1(2,1))+nd3(1,1)*(nd1(2,1)-nd2(2,1)))/2);