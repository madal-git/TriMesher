function NODE=nddata(CELL,M,X,Y,N_L,N_H,N,N_I,N_I_N,h,dx,dy,FM,K)
% NODE=nddata(CELL,M,X,Y,N_L,N_H,N,N_I,N_I_N,h,dx,dy,FM,K) generate the node data structure
% and fill it with basic geometrical information.
% CELL is cell structure data.
% M is the total number of cells.
% X and Y is vector that stores the x and y coordinates of all nodes
% N_L is the number of nodes on the length edge
% N_H is the number of nodes on the height edge
% N is the total number of nodes
% N_I is the total number of interior nodes
% N_I_N is the number of immersed boundary nodes
% h is the node spacing on the outer rectangular boundaries. (For IRT mesh only, h=0 for unstructured mesh)
% dx is the node spacing on the horizontal outer rectangular boundaries. (For random mesh only, h=0)
% dy is the node spacing on the vertical outer rectangular boundaries. (For random mesh only, h=0)
% FM is the flag for mesh type. FM=0----IRT mesh; FM=1----Random mesh
% K is the dimension in the node data structure.

%%%%%%%%%%%%%Structure of each element ND in NODE
%%%%%%%%%%% ND{1}~ND{11} are filled in the current function, where ND{8}~ND{11} for
% corner nodes are filled in bc.m since their values depends on the acutal
% type of boundary condition. In addition, the real more detailed value in
% ND{2} is also filled in bc.m since the ID is BC-dependent.

% ND{1}  The order number of current node, succeed from X and Y
% ND{2}  The identifier of boundary node, 1---On boundary; 0---Interior; -1---Immersed boundary
% ND{3}  Column vector of coordinate of current node

% ND{4}  The total number of triangles connected to current node
% ND{5}  The vector of order number of each connected triangle
% ND{6}  The vector of distance from the centroid each connected triangle
% ND{7}  The sum of reverse distance to each star triangle

% ND{8}  --- The total number of triangles connected to current boundary
% node that could be periodic (Empty for interior nodes)
% ND{9}  --- The vector of order number of each connected triangle for
% periodic boundary node (Empty for interior nodes)
% ND{10} --- The vector of distance from the centroid of each connected
% triangle to periodic boundary node (Empty for interior nodes)
% ND{11} --- The sum of reverse distance from each star triangle to
% periodic boundary node (Empty for interior nodes)

%%%%%%%%%% ND{12}~ND{15} are filled when face data strucure is available
% ND{12}  The total number of faces connected to current node
% ND{13}  The vector of order number of each connected face
%%%%% Nodal star structure for boundary nodes that is potentially periodic %%%%%
% ND{14}  --- The total number of faces connected to current periodic node(Empty for interior nodes or under the third scenario for corner nodes in ND{16})
% ND{15}  --- The vector of order number of each connected face for periodic  boundary node (Empty for interior nodes or under the third scenario for corner nodes in ND{16})

%%%%%%%%%% ND{16}~ND{18} are filled in function stencil.m in the solver
%%%%%%%%%%%%%% Boundary normal stencil for IRT mesh ONLY %%%%%%%%%%%%%%
% ND{16} --- The column vector of order # of the stencil for plat boundaries(Empty for interior nodes)
% [Exterior periodic node #;On-boundary periodic node #;Interior neighbor node #; Interior interior neighbor node #;..]
% ND{17} --- The column vector of order # of the stencil for corner nodes along the diagonal direction (Empty for non-corner nodes).
% [near node #; further node #;..].


% ND{18} --- The clockwise upstream and downstream boundary neighbor node of current boundary node.
% A column vector with first row for upstream node and second row the downstream node(Empty for interior nodes)

%%%%%%%%%% ND{19}~ND{20} are filled in function bc.m in the solver
% ND{19} --- The vector of upwind and bounce-back pdf order # of current node on boundary for V1
% ND{20} --- The vector of upwind and bounce-back pdf order # of current node on boundary for V2

%%%%%%%%%% ND{21} is filled in function macro_bc.m in the solver
% ND{21} --- Dirichlet boundary conditions [density;x velocity;y velocity;temperature]

if K<25
    error('More dimensions in the node structure is required. Please redefine its dimension no less than 25');
end

NODE=cell(K,1);
NODE={NODE,NODE};

X1=min(X);
X2=max(X);
Y1=min(Y);
Y2=max(Y);

P=CELL{1};
e=ceil(10*sqrt(M*P{6})); %%%% Calculate the reference size

if FM==0
    %%%% Find the star structure of cells attached to current node and node
    %%%% information, ND{1}~ND{7}
    %% Method 1--too slow
    % for l=1:N
    %     ND=NODE{l};
    %     ND{1,1}=l;
    %     a=0;
    %     NCM=0;
    %     NCD=0;
    %     DC=0;
    %     for s=1:M
    %         P=CELL{s};
    %         for g=1:3
    %             if P{6+g}==l
    %                 a=a+1;
    %                 NCM(a,1)=P{9+g};
    %                 NCD(a,l)=s;
    %                 Centroid=P{5};
    %                 DC(a,l)=sqrt((X(l)-Centroid(1,1))^2+(Y(l)-Centroid(2,1))^2);
    %                 break;
    %             end
    %         end
    %     end
    %     for g=1:a-1;
    %         if NCM(g,1)~=NCM(g+1,1)
    %             error('The boundary markers conflict for current star triangles!');
    %         end
    %     end
    %     ND{2,1}=NCM(1,1);
    %     ND{3,1}=[X(l);Y(l)];
    %     ND{4,1}=a;
    %     ND{5,1}=NCD(1:a,l);  % 1:a is necessary
    %     ND{6,1}=DC(1:a,l);   % 1:a is necessary
    %     if ND{4}~=length(ND{5}) || ND{4}~=length(ND{6})
    %         error('The length of vector of order # of connected cell is wrong!');
    %     end
    %     D_sum=0;
    %     for g=1:ND{4,1};
    %         D_sum=D_sum+1/DC(g,l);
    %     end
    %     ND{7,1}=D_sum;
    %     NODE{l}=ND;
    % end
    %% Method 2
    N_cell_star=zeros(1,N);
    Order_cell_star=cell(1,N);
    Boundary_ID=cell(1,N);
    for l=1:N
        Order_cell_star{l}=0;
        Boundary_ID{l}=0;
    end
    for r=1:M
        P=CELL{r};
        for g=1:3
            %%% Order number vector of cells
            N_cell_star(P{6+g})=N_cell_star(P{6+g})+1;
            Order_vector=Order_cell_star{P{6+g}};
            Order_vector(1,end+1)=P{1};
            Order_cell_star{P{6+g}}=Order_vector;
            %%% Boundary ID
            B_ID=Boundary_ID{P{6+g}};
            B_ID(end+1)=P{9+g};
            Boundary_ID{P{6+g}}=B_ID;
        end
        
    end
    for l=1:N
        %%% Order number vector of cells
        Order_vector=Order_cell_star{l};
        Order_vector=Order_vector(2:end); % Eliminate the first 0 element
        Order_cell_star{l}=Order_vector;
        %%% Boundary ID
        B_ID=Boundary_ID{l};
        B_ID=B_ID(2:end); % Eliminate the first 0 element
        B_ID=unique(B_ID);
        if length(B_ID)~=1
            error('The same node has multiple boundary identifiers!');
        end
        Boundary_ID{l}=B_ID;
    end
    %%% Create the node data structure
    for l=1:N
        ND=cell(K,1);
        ND{1,1}=l;
        a=N_cell_star(l);
        NCD=Order_cell_star{l};
        for s=1:a
            P=CELL{NCD(s)};
            Centroid=P{5};
            DC(s,l)=sqrt((X(l)-Centroid(1,1))^2+(Y(l)-Centroid(2,1))^2);
        end
        ND{2,1}=Boundary_ID{l};
        ND{3,1}=[X(l);Y(l)];
        ND{4,1}=a;
        ND{5,1}=NCD(1:a);
        ND{6,1}=DC(1:a,l);
        if ND{4}~=length(ND{5}) || ND{4}~=length(ND{6})
            error('The length of vector of order # of connected cell is wrong!');
        end
        D_sum=0;
        for g=1:ND{4,1};
            D_sum=D_sum+1/DC(g,l);
        end
        ND{7,1}=D_sum;
        NODE{l}=ND;
    end
    %% check legitimacy of node data
    n2=0;
    n4=0;
    n8=0;
    n_i=0;
    n_b=0;
    for l=1:N
        ND=NODE{l};
        if ND{4}==2;
            n2=n2+1;
            n_b=n_b+1;
            if ND{2}~=1
                error('!');
            end
        end
        if ND{4}==4;
            n4=n4+1;
            if ND{2}==1
                n_b=n_b+1;
            elseif ND{2}==0
                n_i=n_i+1;
            else
                error('?');
            end
        end
        if ND{4}==8;
            n8=n8+1;
            n_i=n_i+1;
            if ND{2}~=0
                error('!');
            end
        end
    end
    if n2+n4+n8~=N || n_b+n_i~=N
        error('/');
    end
    %% %%%%%% Find the periodic star structure of the nodes on boudaries in order to apply
    % %%%%%% periodic boundary conditions
    % %%%%% For any mesh as long as boundary nodes are lined up
    % % ND{8}  --- The total number of triangles connected to current periodic node
    % % ND{9}  --- The vector of order number of each connected triangle for
    % % periodic boundary node (Empty for interior nodes)
    % % ND{10} --- The vector of distance from the centroid of each connected
    % % triangle to periodic boundary node (Empty for interior nodes)
    % % ND{11} --- The sum of reverse distance from each star triangle to
    % % periodic boundary node (Empty for interior nodes)
    for l=1:N
        if l<=N_I
            ; % Empty for interior nodes
        elseif l==N_I+N_L-1 || l==N_I+N_L-1+N_H-1 || l==N_I+N_L-1+N_H-2+N_L || l==N % Give a temperay value,to be updated in bc.m
            % Each of corner node has eight connected triangles
            NP1=NODE{N_I+N_L-1}; % Top-Right
            NP2=NODE{N_I+N_L-1+N_H-1}; % Bottom-Right
            NP3=NODE{N_I+N_L-1+N_H-2+N_L}; % Bottom-Left
            NP4=NODE{N}; % Top-Left
            ND=NODE{l};
            %%%% Total number of connected cells
            ND{8,1}=NP1{4,1}+NP2{4,1}+NP3{4,1}+NP4{4,1};
            if ND{8,1}~=8
                error('The length of vector of order # of connected cell for corner node is wrong!');
            end
            %%%% Vector of order # of connected cells
            NC=[NP1{5,1},NP2{5,1},NP3{5,1},NP4{5,1}];
            ND{9,1}=NC;
            if ND{8,1}~=length(ND{9,1})
                error('The length of vector of order # of connected cell is wrong!');
            end
            
            % Keep ND{10,1} empty
            % Keep ND{11,1} empty
            %%%% load back to NODE
            NODE{l}=ND;
        else
            ND=NODE{l};
            %%%%% Find the mirror node on the periodic boundary
            if X(l)==X1
                X_target=X2;
                Y_target=Y(l);
            elseif X(l)==X2
                X_target=X1;
                Y_target=Y(l);
            elseif Y(l)==Y1
                X_target=X(l);
                Y_target=Y2;
            elseif Y(l)==Y2
                X_target=X(l);
                Y_target=Y1;
            else
                error('The current node is not on the outer flat boundary (except for the corner nodes) !')
            end
            for t=1:N
                if single(e+X(t))==single(e+X_target) && single(e+Y(t))==single(e+Y_target)
                    NP=NODE{t};
                    break;
                end
            end
            if t==N
                error('mirror node on periodic boundary is not found!');
            end
            %%%% Total number of connected cells
            ND{8,1}=ND{4,1}+NP{4,1};
            %%%% Vector of order # of connected cells
            NC=ND{5,1};
            L=length(NP{5,1});
            NC(end+1:end+L)=NP{5,1};
            ND{9,1}=NC;
            if ND{8,1}~=length(ND{9,1})
                error('The length of vector of order # of connected cell is wrong!');
            end
            %%%% Vector of distance
            NDC=ND{6,1};
            L=length(NP{6,1});
            NDC(end+1:end+L)=NP{6,1};
            ND{10,1}=NDC;
            if ND{8,1}~=length(ND{10,1})
                error('The length of vector of distance is wrong!');
            end
            %%%% Sum of reverse distance
            D_sum=0;
            for g=1:ND{8,1};
                D_sum=D_sum+1/NDC(g);
            end;
            ND{11,1}=D_sum;
            %%%% load back to NODE
            NODE{l}=ND;
        end
    end
elseif FM==1
    disp('This part is moved to PostTri_nddata.m!');
else
    error('Wrong flag for mesh type!');
end