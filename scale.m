function [X,Y,X_B,Y_B,X_I_B,Y_I_B,X_I_AB,Y_I_AB,CELL,h,X1,Y1,X2,Y2]=scale(X,Y,X_B,Y_B,X_I_B,Y_I_B,X_I_AB,Y_I_AB,CELL,M,h,X1,Y1,X2,Y2,S)
% function [X,Y,X_B,Y_B,X_I_B,Y_I_B,X_I_AB,Y_I_AB,CELL,h,X1,Y1,X2,Y2]=scale(X,Y,X_B,Y_B,X_I_B,Y_I_B,X_I_AB,Y_I_AB,CELL,M,h,X1,Y1,X2,Y2,S)
% scales the dimension of all geometric information into real space.

% Scale triangles.
for i=1:M
    P=CELL{i};
    
    CT1=P{13};
    CT1(1,1)=(CT1(1,1)-X1)*S;   % Coordinate of node 1
    CT1(2,1)=(CT1(2,1)-Y1)*S;
    P{13}=CT1;
    
    CT2=P{14};
    CT2(1,1)=(CT2(1,1)-X1)*S;   % Coordinate of node 2
    CT2(2,1)=(CT2(2,1)-Y1)*S;
    P{14}=CT2;

    CT3=P{15};
    CT3(1,1)=(CT3(1,1)-X1)*S;   % Coordinate of node 2
    CT3(2,1)=(CT3(2,1)-Y1)*S;
    P{15}=CT3;
    
    x1=[CT1(1,1),CT2(1,1)];
    x2=[CT2(1,1),CT3(1,1)];
    x3=[CT3(1,1),CT1(1,1)];
    y1=[CT1(2,1),CT2(2,1)];
    y2=[CT2(2,1),CT3(2,1)];
    y3=[CT3(2,1),CT1(2,1)];    

    P{22}=x1;   % x coordinate of face 1
    P{23}=y1;   % y coordinate of face 1

    P{24}=x2;   % x coordinate of face 2
    P{25}=y2;   % y coordinate of face 2

    P{26}=x3;   % x coordinate of face 3
    P{27}=y3;   % y coordinate of face 3

    Centroid=P{5};
    Centroid(1,1)=(CT1(1,1)+CT2(1,1)+CT3(1,1))/3; % Centroid coordinate of current Cell
    Centroid(2,1)=(CT1(2,1)+CT2(2,1)+CT3(2,1))/3;
    P{5}=Centroid;

    Mid1=P{28};
    Mid1(1,1)=(CT1(1,1)+CT2(1,1))/2; % Mid point coordinate of face 1
    Mid1(2,1)=(CT1(2,1)+CT2(2,1))/2;
    P{28}=Mid1;

    Mid2=P{29};
    Mid2(1,1)=(CT2(1,1)+CT3(1,1))/2;    % Mid point coordinate of face 2
    Mid2(2,1)=(CT2(2,1)+CT3(2,1))/2;
    P{29}=Mid2;

    Mid3=P{30};
    Mid3(1,1)=(CT3(1,1)+CT1(1,1))/2;    % Mid point coordinate of face 3
    Mid3(2,1)=(CT3(2,1)+CT1(2,1))/2;
    P{30}=Mid3;
    
    P{6}=abs((CT1(1,1)*(CT2(2,1)-CT3(2,1))+CT2(1,1)*(CT3(2,1)-CT1(2,1))+CT3(1,1)*(CT1(2,1)-CT2(2,1)))/2);   % area of cell
    
    CELL{i}=P;
end;

% Scale all nodes
X=(X-X1)*S;
Y=(Y-Y1)*S;

% Scale boundary nodes
X_B=(X_B-X1)*S;
Y_B=(Y_B-Y1)*S;

% Scale interior nodes that are aligned with boundary nodes
X_I_B=(X_I_B-X1)*S;
Y_I_B=(Y_I_B-Y1)*S;

% Scale interior nodes that are not aligned with boundary nodes
X_I_AB=(X_I_AB-X1)*S;
Y_I_AB=(Y_I_AB-Y1)*S;

% Scale spacing size
h=h*S;

% Scale boundaries
X2=(X2-X1)*S;
Y2=(Y2-Y1)*S;
X1=(X1-X1)*S;
Y1=(Y1-Y1)*S;