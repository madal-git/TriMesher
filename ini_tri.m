function [C,X1,X2,Y1,Y2]=ini_tri(L,RatioLH,K)
% function C=Ini_tri(L,RatioLH,K) defines the initial triangle
% and the cell data structure for following growth.
% L is defined in main function.
% RatioLH is defined in main function.
% K is defined in the main function.
% C is the cell structure containing all mesh information.
% x=X1 is the left face of the rectangular.
% x=X2 is the right face of the rectangular.
% y=Y1 is the bottom face of the rectangular.
% y=Y2 is the top face of the rectangular.

if K<40
    error('More dimensions in the cell structure is required. Please redefine its dimension no less than 40');
end

T1=[0;0];  % Coordinate of point 1 of initial triangle
T2=[5;5];  % Coordinate of point 1 of initial triangle
T3=[10;0]; % Coordinate of point 1 of initial triangle
x1=[T1(1),T2(1)];  % x coordinate of edge 1
x2=[T2(1),T3(1)];  % x coordinate of edge 2
x3=[T3(1),T1(1)];  % x coordinate of edge 3
y1=[T1(2),T2(2)];  % y coordinate of edge 1
y2=[T2(2),T3(2)];  % y coordinate of edge 2
y3=[T3(2),T1(2)];  % y coordinate of edge 3

C=cell(K,1);

C{1}=1;  % Cell Number

C{2}=0;  % The order number of neighbor triangle on face 1 (0:No neighbor; nonezero:has neighbor)
C{3}=0;  % The order number of neighbor triangle on face 2 (0:No neighbor; nonezero:has neighbor)
C{4}=0;  % The order number of neighbor triangle on face 3 (0:No neighbor; nonezero:has neighbor)

C{5}=[(T1(1)+T2(1)+T3(1))/3;(T1(2)+T2(2)+T3(2))/3]; % Centroid coordinate of Cell

C{6}=abs((T1(1)*(T2(2)-T3(2))+T2(1)*(T3(2)-T1(2))+T3(1)*(T1(2)-T2(2)))/2);   % area of cell

C{7}=1;    % Number of node 1
C{8}=2;    % Number of node 2
C{9}=3;    % Number of node 3

C{10}=1;   % Region identifier of node 1 (1---Wall Boundary; 0---Interior)
C{11}=1;   % Region identifier of node 2 (1---Wall Boundary; 0---Interior)
C{12}=1;   % Region identifier of node 3 (1---Wall Boundary; 0---Interior)

C{13}=[T1(1);T1(2)];   % Coordinate of node 1
C{14}=[T2(1);T2(2)];   % Coordinate of node 2
C{15}=[T3(1);T3(2)];   % Coordinate of node 3

C{16}=1;   % Order Number of face 1
C{17}=2;   % Order Number of face 2
C{18}=3;   % Order Number of face 3

C{19}=1;   % Region identifier of face 1 (1---Wall Boundary; 0---Interior)
C{20}=1;   % Region identifier of face 2 (1---Wall Boundary; 0---Interior)
C{21}=1;   % Region identifier of face 3 (1---Wall Boundary; 0---Interior)

C{22}=x1;   % x coordinate of face 1
C{23}=y1;   % y coordinate of face 1
C{24}=x2;   % x coordinate of face 2
C{25}=y2;   % y coordinate of face 2
C{26}=x3;   % x coordinate of face 3
C{27}=y3;   % y coordinate of face 3

C{28}=[(T1(1)+T2(1))/2;(T1(2)+T2(2))/2];    % Mid point coordinate of face 1
C{29}=[(T2(1)+T3(1))/2;(T2(2)+T3(2))/2];    % Mid point coordinate of face 2
C{30}=[(T3(1)+T1(1))/2;(T3(2)+T1(2))/2];    % Mid point coordinate of face 3

C{31}=sqrt(50);   % The length of face 1
C{32}=sqrt(50);   % The length of face 2
C{33}=10;         % The length of face 3

C{34}=[-sqrt(2)/2; sqrt(2)/2];   % The unit normal vector of face 1
C{35}=[sqrt(2)/2; sqrt(2)/2];    % The unit normal vector of face 2
C{36}=[0;1];                     % The unit normal vector of face 3

C={C,C};

X1=T2(1)-L/2;
X2=T2(1)+L/2;
Y1=T2(2)*(1/3)-0.2*(L-2)-L/2/RatioLH;
Y2=T2(2)*(1/3)-0.2*(L-2)+L/2/RatioLH;