function [C,M]=IRTtriangulation(N_H,N_L,C,X,Y,N,X_I_AB,Y_I_AB,h,X1,Y1,X2,Y2,K)
% [C,M]=IRTtriangulation(N_H,N_L,C,X,Y,N,X_I_AB,h,Y_I_AB) generate
% isosceles right triangulation.
% N_H is defined in the main function.
% N_L is defined in the function B_node.
% C is input and ouput in the same time.
% X, Y and N are defined in the function node_combine.
% X_I_AB and Y_I_AB are defined in the function I_nodeAB.
% h is defined in function B_node.
% X1,Y1,X2,Y2 are defined in the function ini_tri.
% M is the total number of triangles plus the initial big triangle.
% K is total number of data input for each cell data.

if K<40
    error('More dimensions in the cell structure is required. Please redefine its dimension no less than 40');
end

% Phase pointer
OR(1,1:4)=[-1,1,1,1];
OR(2,1:4)=[1,1,1,-1];
OR(3,1:4)=[1,-1,-1,-1];
OR(4,1:4)=[-1,-1,-1,1];
% Triangulation starts
M=1;
for i=1:N_L-1
    for j=1:N_H-1
        for g=1:4
            M=M+1;
            CC=cell(K,1);
            CC{1}=M;  % Cell Number

            CC{10}=0;   % Region identifier of node 1 (1---Wall Boundary; 0---Interior)
            CC{13}=[X_I_AB(i,j);Y_I_AB(i,j)];   % Coordinate of node 1
            C_T_1=CC{13};
            for n=1:N;
                if single(C_T_1(1,1))==single(X(n)) && single(C_T_1(2,1))==single(Y(n));
                    break;
                end
            end
            CC{7}=n;   % Number of node 1

            CC{11}=0;    % Region identifier of node 2 (1---Wall Boundary; 0---Interior)
            CC{14}=[X_I_AB(i,j)+OR(g,1)*h/2;Y_I_AB(i,j)+OR(g,2)*h/2];   % Coordinate of node 2
            C_T_2=CC{14};
            for n=1:N
                if single(C_T_2(1,1))==single(X(n)) && single(C_T_2(2,1))==single(Y(n));
                    break;
                end
            end
            CC{8}=n;    % Number of node 2

            CC{12}=0;    % Region identifier of node 3 (1---Wall Boundary; 0---Interior)
            CC{15}=[X_I_AB(i,j)+OR(g,3)*h/2;Y_I_AB(i,j)+OR(g,4)*h/2];   % Coordinate of node 3
            C_T_3=CC{15};
            for n=1:N
                if single(C_T_3(1,1))==single(X(n)) && single(C_T_3(2,1))==single(Y(n));
                    break;
                end
            end
            CC{9}=n;    % Number of node 3
                        
            x1=[C_T_1(1,1),C_T_2(1,1)];
            x2=[C_T_2(1,1),C_T_3(1,1)];
            x3=[C_T_3(1,1),C_T_1(1,1)];
            y1=[C_T_1(2,1),C_T_2(2,1)];
            y2=[C_T_2(2,1),C_T_3(2,1)];
            y3=[C_T_3(2,1),C_T_1(2,1)];
            CC{22}=x1;   % x coordinate of face 1
            CC{23}=y1;   % y coordinate of face 1

            CC{24}=x2;   % x coordinate of face 2
            CC{25}=y2;   % y coordinate of face 2

            CC{26}=x3;   % x coordinate of face 3
            CC{27}=y3;   % y coordinate of face 3

            CC{5}=[(C_T_1(1,1)+C_T_2(1,1)+C_T_3(1,1))/3;(C_T_1(2,1)+C_T_2(2,1)+C_T_3(2,1))/3]; % Centroid coordinate of current Cell

            CC{28}=[(C_T_1(1,1)+C_T_2(1,1))/2;(C_T_1(2,1)+C_T_2(2,1))/2];    % Mid point coordinate of face 1

            CC{29}=[(C_T_2(1,1)+C_T_3(1,1))/2;(C_T_2(2,1)+C_T_3(2,1))/2];    % Mid point coordinate of face 2

            CC{30}=[(C_T_3(1,1)+C_T_1(1,1))/2;(C_T_3(2,1)+C_T_1(2,1))/2];    % Mid point coordinate of face 3

            CC{6}=area_tri(C_T_1,C_T_2,C_T_3);   % area of cell
                
            C{M}=CC;
%             Message=['Triangle ' , num2str(M-1) , ' is created'];
%             disp(Message);
        end
    end
end
%%%%% Define the new boundary identifier
for p=1:M
    T=C{p};
    CT1=T{13};
    CT2=T{14};
    CT3=T{15};
    if single(CT1(1,1))==single(X1) || single(CT1(1,1))==single(X2) || single(CT1(2,1))==single(Y1) || single(CT1(2,1))==single(Y2)
        T{10}=1;
    end
    
    if single(CT2(1,1))==single(X1) || single(CT2(1,1))==single(X2) || single(CT2(2,1))==single(Y1) || single(CT2(2,1))==single(Y2)
        T{11}=1;
    end
    
    if single(CT3(1,1))==single(X1) || single(CT3(1,1))==single(X2) || single(CT3(2,1))==single(Y1) || single(CT3(2,1))==single(Y2)
        T{12}=1;
    end   
    C{p}=T;
end
% Get rid of initial triangle
M=M-1;