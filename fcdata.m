function [FACE,O]=fcdata(CELL,NODE,X,Y,X_B,Y_B,N_B,X_I_B,Y_I_B,N_I_B,X_I_AB,Y_I_AB,N_I_AB,N_L,N_H,X1,X2,Y1,Y2,h,FM,K)
% [FACE,O]=fcdata(CELL,NODE,X1,X2,Y1,Y2,FM) generates the face data structure with
% basic geometrical information.
% FACE is the face data structure for all faces.
% O is the total number of faces.
% CELL is the cell data structure for all triangles.
% NODE is the node data structure for all nodes.
% X1, X2, Y1 and Y2 are bounds of mesh domain
% FM is flag for mesh type. FM=0----IRT mesh; FM=1----Random mesh
% K is the dimension in the face data structure.

% FACE is the resulting data structure. the information stored for each
% face is listed as follows.

% FC=FACE{r}

%%%%%%%% FC{1}~{11} will be filled in the current function,where FC{2},
%%%%%%%% FC{10}~{11} will be updated in bc.m
% FC{1}   Current face order #
% FC{2}   Boundary identifier of current face(0-Interior;1-Periodic;20-Velocity Inlet;21-Pressure Inlet;30-Velocity Outlet;31-Pressure Outlet;4-S Wall;5-M wall;6-Zero gradient. Negative for Immersed boundaries)
% FC{3}   The length of current face
% FC{4}   The unit normal vector of current face (right-hand-rule)
% FC{5}   The coordinate of first end node, column vector [x;y]
% FC{6}   The coordinate of first end node, column vector [x;y]
% FC{7}   The coordinate of mid point, column vector [x;y]

% FC{8}   The order # of 1st end node
% FC{9}   The order # of 2nd end node
% FC{10}  Boundary identifier of 1st end node
% FC{11}  Boundary identifier of 2nd end node

%%%%%%%% FC{12}~{15} will be filled in IRTngb.m
% FC{12}  The order # of 1st neighbor cell (0:No neighbor; nonezero:has neighbor)
% FC{13}  The order # of 2nd neighbor cell (0:No neighbor; nonezero:has neighbor)

%%%%%%%%%%%%%%%%%%%%%%%%%% For Periodic boundaries ONLY %%%%%%%%%%%%%%%%%%%%%%%%%%%
% FC{14}  The order # of 1st neighbor cell (nonezero and FC{2} has to be 1)
% FC{15}  The order # of 2nd neighbor cell (nonezero and FC{2} has to be 1)
%%%%%%%%%%%%%%%%%%%%%%%%%% For Periodic boundaries ONLY %%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% FC{16}~{23} will be filled in stencil.m
% FC{16}   The vector of triangle order # of 1st upwind scheme of current face for V1
% FC{17}   The vector of triangle order # of 1st upwind scheme of current face for V2
% FC{18}   The matrix of triangle order # of 2nd upwind scheme of current face for V1
% FC{19}   The matrix of triangle order # of 2nd upwind scheme of current face for V2
%% for FC{18} and FC{19},each column is
% [1 downwind cell #;
%  2 upwind cell #;
%  3 further upwind cell #;
%  4 distance from downwind cell to upwind cell;
%  5 distance from upwind cell to further upwind cell;
%  6 spatial ratio for current edge,R1=(downwind cell size + upwind cell size)/upwind cell size;
%  7 spatial ratio for extrapolation, R2=distance from downwind cell to upwind cell/ distance from upwind cell to further upwind cell
%  8 the left boundary node # on edge for extrapolation;
%  9 the right boundary node # on edge for extrapolation;
%  10 spatial ratio for on edge,R3=distance from left node to intercept node/distance from left node to right node]

%%%%%%%%%%%%%%%%%%%%%%%%%% For Periodic boundaries ONLY %%%%%%%%%%%%%%%%%%%%%%%%%%%
% FC{20}  The vector of triangle order # of 1st upwind of current face on boundary for V1
% FC{21}  The vector of triangle order # of 1st upwind of current face on boundary for V2

% FC{22}  The matrix of triangle order # of 2nd upwind of current face on boundary for V1
% FC{23}  The matrix of triangle order # of 2nd upwind of current face on boundary for V2
%%%%%%%%%%%%%%%%%%%%%%%%%% For Periodic boundaries ONLY %%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% FC{24}~{25} will be filled in pdf_bounceback.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%% flux bounceback %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FC{24}  The vector of pdf order # of upwind and bounceback of current face on boundary for V1
% FC{25}  The vector of pdf order # of upwind and bounceback of current face on boundary for V2
%%%%%%%%%%%%%%%%%%%%%%%%%%%% flux bounceback %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if K<55
    error('More dimensions in the face data structure is required. Please redefine its dimension no less than 55');
end

O=0; % counter for faces
FACE=cell(K,1);
FACE={FACE,FACE};
M=length(CELL);
N=length(NODE);
P=CELL{1};
e=ceil(10*sqrt(M*P{6})); %%%% Calculate the reference size

if FM==0
    %% Phase pointer
    OR(1,1:2)=[-1, 1];
    OR(2,1:2)=[ 1, 1];
    OR(3,1:2)=[ 1,-1];
    OR(4,1:2)=[-1,-1];
    %% Generate the faces that share the nodes which are not aligned with the boundary nodes
    for i=1:N_L-1
        for j=1:N_H-1
            % The coordinates of first node
            C_nd1=[X_I_AB(i,j);Y_I_AB(i,j)];   % Coordinate of node 1
            % find the order number of the first node
            for s=1:N
                if single(e+C_nd1(1,1))==single(e+X(s)) && single(e+C_nd1(2,1))==single(e+Y(s));
                    break;
                end
            end
            if s==N
                if single(e+C_nd1(1,1))~=single(e+X(s)) || single(e+C_nd1(2,1))~=single(e+Y(s));
                    error('The order number for the first node that is not aligned with the boundary node is not found!');
                end
            end
            ND1=NODE{s};
            for g=1:4
                O=O+1;
                FC=cell(K,1);
                FC{1}=O;  % Cell Number
                FC{2}=0;  % boundary identifier
                
                FC{8}=ND1{1}; % The order number of first node
                FC{10}=ND1{2}; % The order number of first node
                
                C_nd2=[X_I_AB(i,j)+OR(g,1)*h/2;Y_I_AB(i,j)+OR(g,2)*h/2];   % Coordinate of node 2
                for n=1:N
                    if single(e+C_nd2(1,1))==single(e+X(n)) && single(e+C_nd2(2,1))==single(e+Y(n));
                        break;
                    end
                end
                if n==N
                    if single(e+C_nd2(1,1))~=single(e+X(n)) || single(e+C_nd2(2,1))~=single(e+Y(n));
                        error('The order number for the second node that is aligned with the boundary node is not found!');
                    end
                end
                ND2=NODE{n};
                FC{9}=ND2{1}; % The order number of first node
                FC{11}=ND2{2}; % The order number of first node
                
                [Le,Ne]=en(C_nd1,C_nd2);
                FC{3}=Le;
                FC{4}=Ne;
                FC{5}=C_nd1;
                FC{6}=C_nd2;
                FC{7}=(C_nd1+C_nd2)/2;
                
                FACE{O}=FC;
            end
        end
    end
    %% Generate the vertical faces that share the nodes which are aligned with the boundary nodes
    for i=1:N_L-2
        for j=1:N_H-1
            O=O+1;
            % Determine ND1 and ND2
            if j==N_H-1
                C_nd1=[X_I_B(i,N_H-2);Y_I_B(i,N_H-2)];   % This will cause round-off difference from ND{3}
                C_nd2=[X_I_B(i,N_H-2);Y_I_B(i,N_H-2)+h]; % This will cause round-off difference from ND{3}
            else
                C_nd1=[X_I_B(i,j);Y_I_B(i,j)-h]; % This will cause round-off difference from ND{3}
                C_nd2=[X_I_B(i,j);Y_I_B(i,j)]; % This will cause round-off difference from ND{3}
            end
            for s=1:N
                if single(e+C_nd1(1,1))==single(e+X(s)) && single(e+C_nd1(2,1))==single(e+Y(s));
                    break;
                end
            end
            if s==N
                if single(e+C_nd1(1,1))~=single(e+X(s)) || single(e+C_nd1(2,1))~=single(e+Y(s));
                    error('The order number for the first node that is not aligned with the boundary node is not found!');
                end
            end
            ND1=NODE{s};
            for n=1:N
                if single(e+C_nd2(1,1))==single(e+X(n)) && single(e+C_nd2(2,1))==single(e+Y(n));
                    break;
                end
            end
            if n==N
                if single(e+C_nd2(1,1))~=single(e+X(n)) || single(e+C_nd2(2,1))~=single(e+Y(n));
                    error('The order number for the second node that is aligned with the boundary node is not found!');
                end
            end
            ND2=NODE{n};
            % Fill the data
            FC=cell(K,1);
            FC{1}=O;  % Cell Number
            FC{2}=0;  % boundary identifier
            
            FC{8}=ND1{1}; % The order number of first node
            FC{10}=ND1{2}; % The order number of first node
            FC{9}=ND2{1}; % The order number of first node
            FC{11}=ND2{2}; % The order number of first node
            
            [Le,Ne]=en(C_nd1,C_nd2);
            FC{3}=Le;
            FC{4}=Ne;
            FC{5}=C_nd1;
            FC{6}=C_nd2;
            FC{7}=(C_nd1+C_nd2)/2;
            
            FACE{O}=FC;
        end
    end
    %% Generate the horizontal faces that share the nodes which are aligned with the boundary nodes
    for j=1:N_H-2
        for i=1:N_L-1
            O=O+1;
            % Determine ND1 and ND2
            if i==N_L-1
                C_nd1=[X_I_B(N_L-2,j);Y_I_B(N_L-2,j)]; % This will cause round-off difference from ND{3}
                C_nd2=[X_I_B(N_L-2,j)+h;Y_I_B(N_L-2,j)]; % This will cause round-off difference from ND{3}
            else
                C_nd1=[X_I_B(i,j)-h;Y_I_B(i,j)]; % This will cause round-off difference from ND{3}
                C_nd2=[X_I_B(i,j);Y_I_B(i,j)]; % This will cause round-off difference from ND{3}
            end
            for s=1:N
                if single(e+C_nd1(1,1))==single(e+X(s)) && single(e+C_nd1(2,1))==single(e+Y(s));
                    break;
                end
            end
            if s==N
                if single(e+C_nd1(1,1))~=single(e+X(s)) || single(e+C_nd1(2,1))~=single(e+Y(s));
                    error('The order number for the first node that is not aligned with the boundary node is not found!');
                end
            end
            ND1=NODE{s};
            for n=1:N
                if single(e+C_nd2(1,1))==single(e+X(n)) && single(e+C_nd2(2,1))==single(e+Y(n));
                    break;
                end
            end
            if n==N
                if single(e+C_nd2(1,1))~=single(e+X(n)) || single(e+C_nd2(2,1))~=single(e+Y(n));
                    error('The order number for the second node that is aligned with the boundary node is not found!');
                end
            end
            ND2=NODE{n};
            % Fill the data
            FC=cell(K,1);
            FC{1}=O;  % Cell Number
            FC{2}=0;  % boundary identifier
            
            FC{8}=ND1{1}; % The order number of first node
            FC{10}=ND1{2}; % The order number of first node
            FC{9}=ND2{1}; % The order number of first node
            FC{11}=ND2{2}; % The order number of first node
            
            [Le,Ne]=en(C_nd1,C_nd2);
            FC{3}=Le;
            FC{4}=Ne;
            FC{5}=C_nd1;
            FC{6}=C_nd2;
            FC{7}=(C_nd1+C_nd2)/2;
            
            FACE{O}=FC;
        end
    end
    %% Generate the faces on the boundary
    for i=1:N_B
        O=O+1;
        % Determine ND1 and ND2
        if i==N_B
            C_nd1=[X_B(i);Y_B(i)];
            C_nd2=[X_B(1);Y_B(1)];
        else
            C_nd1=[X_B(i);Y_B(i)];
            C_nd2=[X_B(i+1);Y_B(i+1)];
        end
        for s=1:N
            if single(e+C_nd1(1,1))==single(e+X(s)) && single(e+C_nd1(2,1))==single(e+Y(s));
                break;
            end
        end
        if s==N
            if single(e+C_nd1(1,1))~=single(e+X(s)) || single(e+C_nd1(2,1))~=single(e+Y(s));
                error('The order number for the first node that is not aligned with the boundary node is not found!');
            end
        end
        ND1=NODE{s};
        for n=1:N
            if single(e+C_nd2(1,1))==single(e+X(n)) && single(e+C_nd2(2,1))==single(e+Y(n));
                break;
            end
        end
        if n==N
            if single(e+C_nd2(1,1))~=single(e+X(n)) || single(e+C_nd2(2,1))~=single(e+Y(n));
                error('The order number for the second node that is aligned with the boundary node is not found!');
            end
        end
        ND2=NODE{n};
        % Fill the data
        FC=cell(K,1);
        FC{1}=O;  % Cell Number
        FC{2}=1;  % boundary identifier
        
        FC{8}=ND1{1}; % The order number of first node
        FC{10}=ND1{2}; % The order number of first node
        FC{9}=ND2{1}; % The order number of first node
        FC{11}=ND2{2}; % The order number of first node
        
        [Le,Ne]=en(C_nd1,C_nd2);
        FC{3}=Le;
        FC{4}=Ne;
        FC{5}=C_nd1;
        FC{6}=C_nd2;
        FC{7}=(C_nd1+C_nd2)/2;
        
        FACE{O}=FC;
    end
elseif FM==1
    disp('This part is moved to PostTri_fcdata.m!');
else
    error('The flag for mesh type is incorrect!');
end