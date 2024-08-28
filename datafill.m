function [CELL,NODE,FACE]=datafill(N_I,N_L,N_H,X,Y,CELL,NODE,FACE,FM)
% [CELL,NODE,FACE]=datafill(N_I,N_L,N_H,X,Y,CELL,NODE,FACE,FM) fill information into
% some entries in CELL, NODE and FACE
% CELL is the cell data structure
% NODE is the node data structure
% FACE is the face data structure
% N_I is the total number of interior nodes
% N_L is the number of nodes on the length edge
% N_H is the number of nodes on the height edge
% X and Y is vector that stores the x and y coordinates of all nodes
% FM is the flag for mesh type. FM=0----IRT mesh; FM=1----General mesh

%% For IRT mesh (FM=0), the data to be filled is listed as follows.
%%%%%%%%%% In CELL, C=CELL{m}
% C{16}   Order Number of face 1
% C{17}   Order Number of face 2
% C{18}   Order Number of face 3

% C{19}   Region identifier of face 1
% C{20}   Region identifier of face 2
% C{21}   Region identifier of face 3

% C{31}   The length of face 1
% C{32}   The length of face 2
% C{33}   The length of face 3

% C{34}   The unit normal vector of face 1 (outwards of current cell CV)
% C{35}   The unit normal vector of face 2 (outwards of current cell CV)
% C{36}   The unit normal vector of face 3 (outwards of current cell CV)
%%%%%%%%%% In NODE, ND=NODE{n}
% ND{12}  The total number of faces connected to current node
% ND{13}  The vector of order number of each connected face
%%%%% Nodal star structure for boundary nodes that is potentially periodic %%%%%
% ND{14}  --- The total number of faces connected to current periodic node(Empty for interior nodes or under the third scenario for corner nodes in ND{16})
% ND{15}  --- The vector of order number of each connected face for periodic  

%% For General mesh (FM=1), the data to be filled is listed as follows.
%%%%%%%%%% In CELL, C=CELL{m}



%%%%%%%%%% In NODE, ND=NODE{n}



%% %% %% %% %% %%
X1=min(X);
X2=max(X);
Y1=min(Y);
Y2=max(Y);
M=length(CELL);
N=length(NODE);
O=length(FACE);
P=CELL{1};
e=ceil(10*sqrt(M*P{6})); %%%% Calculate the reference size
if FM==0
    %% Fill C{16}~{21}, C{31}~{36} in the cell data structure
    % Method 1, too slow O(N^2)
%     for r=1:M
%         P=CELL{r};
%         ND_ID1=P{7};
%         ND_ID2=P{8};
%         ND_ID3=P{9};
%         face_found_counter=0;
%         for l=1:O
%             FC=FACE{l};
%             FC_END_ID1=FC{8};
%             FC_END_ID2=FC{9};
%             % Face 1
%             if (FC_END_ID1==ND_ID1 && FC_END_ID2==ND_ID2) || (FC_END_ID1==ND_ID2 && FC_END_ID2==ND_ID1) % Changing
%                 face_found_counter=face_found_counter+1;
%                 % Check legitimacy using face mid point
%                 Mid_cell=P{28}; % changing
%                 Mid_face=FC{7};
%                 if single(e+Mid_cell(1,1))~=single(e+Mid_face(1,1)) || single(e+Mid_cell(2,1))~=single(e+Mid_face(2,1))
%                     error('The found face does not match the current cell!');
%                 end
%                 P{16}=FC{1}; % Changing
%                 P{19}=FC{2}; % Changing
%                 P{31}=FC{3}; % Changing
%                 Centroid=P{5};
%                 Out_n=(Mid_cell-Centroid)'; % The vector that is pointing outwards of current cell control volume
%                 if Out_n*FC{4}'>0
%                     P{34}=FC{4}; % Changing
%                 else
%                     P{34}=-FC{4}; % Changing
%                 end
%             end
%             % Face 2
%             if (FC_END_ID1==ND_ID2 && FC_END_ID2==ND_ID3) || (FC_END_ID1==ND_ID3 && FC_END_ID2==ND_ID2) % Changing
%                 face_found_counter=face_found_counter+1;
%                 % Check legitimacy using face mid point
%                 Mid_cell=P{29}; % changing
%                 Mid_face=FC{7};
%                 if single(e+Mid_cell(1,1))~=single(e+Mid_face(1,1)) || single(e+Mid_cell(2,1))~=single(e+Mid_face(2,1))
%                     error('The found face does not match the current cell!');
%                 end
%                 P{17}=FC{1}; % Changing
%                 P{20}=FC{2}; % Changing
%                 P{32}=FC{3}; % Changing
%                 Centroid=P{5};
%                 Out_n=(Mid_cell-Centroid)'; % The vector that is pointing outwards of current cell control volume
%                 if Out_n*FC{4}'>0
%                     P{35}=FC{4}; % Changing
%                 else
%                     P{35}=-FC{4}; % Changing
%                 end
%             end
%             % Face 3
%             if (FC_END_ID1==ND_ID3 && FC_END_ID2==ND_ID1) || (FC_END_ID1==ND_ID1 && FC_END_ID2==ND_ID3) % Changing
%                 face_found_counter=face_found_counter+1;
%                 % Check legitimacy using face mid point
%                 Mid_cell=P{30}; % changing
%                 Mid_face=FC{7};
%                 if single(e+Mid_cell(1,1))~=single(e+Mid_face(1,1)) || single(e+Mid_cell(2,1))~=single(e+Mid_face(2,1))
%                     error('The found face does not match the current cell!');
%                 end
%                 P{18}=FC{1}; % Changing
%                 P{21}=FC{2}; % Changing
%                 P{33}=FC{3}; % Changing
%                 Centroid=P{5};
%                 Out_n=(Mid_cell-Centroid)'; % The vector that is pointing outwards of current cell control volume
%                 if Out_n*FC{4}'>0
%                     P{36}=FC{4}; % Changing
%                 else
%                     P{36}=-FC{4}; % Changing
%                 end
%             end
%         end
%         if face_found_counter~=3
%             error('The current cell does not has three faces!');
%         end
%         CELL{r}=P;
%     end
    % Method 2, ~O(N)
    for l=1:O
        FC=FACE{l};
        FC_END_ID1=FC{8};
        FC_END_ID2=FC{9};
        ND1=NODE{FC_END_ID1};
        ND2=NODE{FC_END_ID2};
        Cell_star1=ND1{5};
        Cell_star2=ND2{5};
        Common_cell=intersect(Cell_star1,Cell_star2);
        L=length(Common_cell); % L=1 or 2
        if L==0 || L>2
            error('Each face should has one or two cells attached!');
        end
        for i=1:L
            face_found_counter=0;
            P=CELL{Common_cell(i)};
            ND_ID1=P{7};
            ND_ID2=P{8};
            ND_ID3=P{9};
            % Face 1
            if (FC_END_ID1==ND_ID1 && FC_END_ID2==ND_ID2) || (FC_END_ID1==ND_ID2 && FC_END_ID2==ND_ID1) % Changing
                face_found_counter=face_found_counter+1;
                % Check legitimacy using face mid point
                Mid_cell=P{28}; % changing
                Mid_face=FC{7};
                if single(e+Mid_cell(1,1))~=single(e+Mid_face(1,1)) || single(e+Mid_cell(2,1))~=single(e+Mid_face(2,1))
                    error('The found face does not match the current cell!');
                end
                P{16}=FC{1}; % Changing
                P{19}=FC{2}; % Changing
                P{31}=FC{3}; % Changing
                Centroid=P{5};
                Out_n=(Mid_cell-Centroid)'; % The vector that is pointing outwards of current cell control volume
                if Out_n*FC{4}'>0
                    P{34}=FC{4}; % Changing
                else
                    P{34}=-FC{4}; % Changing
                end
            end
            % Face 2
            if (FC_END_ID1==ND_ID2 && FC_END_ID2==ND_ID3) || (FC_END_ID1==ND_ID3 && FC_END_ID2==ND_ID2) % Changing
                face_found_counter=face_found_counter+1;
                % Check legitimacy using face mid point
                Mid_cell=P{29}; % changing
                Mid_face=FC{7};
                if single(e+Mid_cell(1,1))~=single(e+Mid_face(1,1)) || single(e+Mid_cell(2,1))~=single(e+Mid_face(2,1))
                    error('The found face does not match the current cell!');
                end
                P{17}=FC{1}; % Changing
                P{20}=FC{2}; % Changing
                P{32}=FC{3}; % Changing
                Centroid=P{5};
                Out_n=(Mid_cell-Centroid)'; % The vector that is pointing outwards of current cell control volume
                if Out_n*FC{4}'>0
                    P{35}=FC{4}; % Changing
                else
                    P{35}=-FC{4}; % Changing
                end
            end
            % Face 3
            if (FC_END_ID1==ND_ID3 && FC_END_ID2==ND_ID1) || (FC_END_ID1==ND_ID1 && FC_END_ID2==ND_ID3) % Changing
                face_found_counter=face_found_counter+1;
                % Check legitimacy using face mid point
                Mid_cell=P{30}; % changing
                Mid_face=FC{7};
                if single(e+Mid_cell(1,1))~=single(e+Mid_face(1,1)) || single(e+Mid_cell(2,1))~=single(e+Mid_face(2,1))
                    error('The found face does not match the current cell!');
                end
                P{18}=FC{1}; % Changing
                P{21}=FC{2}; % Changing
                P{33}=FC{3}; % Changing
                Centroid=P{5};
                Out_n=(Mid_cell-Centroid)'; % The vector that is pointing outwards of current cell control volume
                if Out_n*FC{4}'>0
                    P{36}=FC{4}; % Changing
                else
                    P{36}=-FC{4}; % Changing
                end
            end
            if face_found_counter~=1
                error('The current face shared multiple times by the same cell, or the current face has no cell sharing it!');
            end
            CELL{Common_cell(i)}=P;
        end
    end
    %% Fill ND{12}~{13} in the node data structure
    N_face_star=zeros(1,N);
    Order_face_star=cell(1,N);
    for l=1:N
        Order_face_star{l}=0;
    end
    for s=1:O
        FC=FACE{s};
        for g=1:2
            %%% Order number vector of cells
            N_face_star(FC{7+g})=N_face_star(FC{7+g})+1;
            Order_vector=Order_face_star{FC{7+g}};
            Order_vector(end+1)=FC{1};
            Order_face_star{FC{7+g}}=Order_vector;
        end
    end
    for l=1:N
        %%% Order number vector of cells
        Order_vector=Order_face_star{l};
        Order_vector=Order_vector(2:end);
        Order_face_star{l}=unique(Order_vector);
        if N_face_star(l)~=length(Order_face_star{l})
            error('The numner of star faces for the current node is incorrect!');
        end
    end
    % Fill the data
    for l=1:N
        ND=NODE{l};
        
        ND{12,1}=N_face_star(l);
        ND{13,1}=Order_face_star{l};

        NODE{l}=ND;
    end
    %% Fill ND{14}~{15} in the node data structure
    for l=1:N
        if l<=N_I
            ; % Empty for interior nodes
        elseif l==N_I+N_L-1 || l==N_I+N_L-1+N_H-1 || l==N_I+N_L-1+N_H-2+N_L || l==N % To be done in bc.m
            ;
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
            %%%% Total number of connected faces
            ND{14,1}=ND{12,1}+NP{12,1};
            %%%% Vector of order # of connected faces
            NC=ND{13,1};
            L=length(NP{13,1});
            NC(end+1:end+L)=NP{13,1};
            ND{15,1}=NC;
            if ND{14,1}~=length(ND{15,1})
                error('The length of vector of order # of connected faces is wrong!');
            end
            %%%% load back to NODE
            NODE{l}=ND;
        end
    end
elseif FM==1
    disp('This part is moved to PostTri_datafill_check.m!');
else
    error('The flag for mesh type is incorrect!');
end