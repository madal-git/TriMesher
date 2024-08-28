function A = tri_area(P1,P2,P3)
% function A = tri_area(P1,P2,P3) calculates the area of a triagnle with
% three vertices P1, P2 and P3
% P1, P2 and P3 are the coordinates of three vertices, which must be column
% vectors and must be 2D
% A is the area

[D1,L1]=size(P1);
[D2,L2]=size(P2);
[D3,L3]=size(P3);

if (L1~=1 || L2~=1) || L3~=1
    error('One of the vertices has multiple or zero coordinates!');
end

if (D1~=2 || D2~=2) || D3~=2
    error('One of the vertices has non-2D coordinates!');
end

A=abs((P1(1)*(P2(2)-P3(2))+P2(1)*(P3(2)-P1(2))+P3(1)*(P1(2)-P2(2)))/2);