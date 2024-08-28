function [FACE,O]=PostTri_fcdata(CELL,NODE,edge,N_I_N,N_I,N_H,N_L,X,Y,X1,X2,Y1,Y2,dx,dy)
% function [FACE,O]=PostTri_fcdata(CELL,NODE,X,Y,N_I_N,N_I,N_L,N_H,X1,X2,Y1,Y2,dx,dy)
% FACE is the FACE data structure for all faces.
% O is the total number of faces.

% CELL is the CELL data structure.
% NODE is the NODE data structure
% edge is the raw data generated from Triangle
% N_I_N is the number of immersed boundary nodes
% N_I is the total number of interior nodes including immersed
% N_H is the number of nodes on the height edge
% N_L is the number of nodes on the length edge
% X and Y is vector that stores the x and y coordinates of all nodes
% X1, X2, Y1 and Y2 are bounds of mesh domain
% dx is the node spacing on the horizontal outer rectangular boundaries.
% dy is the node spacing on the vertical outer rectangular boundaries.

% FACE is the resulting data structure. the information stored for each 
% face is listed as follows.

% FC=FACE{r}

%%%%%%%% FC{1}~{15} will be filled in the current function,where FC{2} and FC{10}~{11} will be updated in bc.m

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

% FC{12}  The order # of 1st neighbor cell (0:No neighbor; nonezero:has neighbor)
% FC{13}  The order # of 2nd neighbor cell (0:No neighbor; nonezero:has neighbor)

%%%%%%%%%%%%%%%%%%%%%%%%%% For Periodic boundaries ONLY %%%%%%%%%%%%%%%%%%%%%%%%%%%
% FC{14}  The order # of 1st neighbor cell (nonezero and FC{2} has to be 1)
% FC{15}  The order # of 2nd neighbor cell (nonezero and FC{2} has to be 1)
%%%%%%%%%%%%%%%%%%%%%%%%%% For Periodic boundaries ONLY %%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%%%%%% FC{16}~{18} are purely mesh info and lattice-independent, because
%%%%%%% of which they should be filled here. But instead, they are filled
%%%%%%% in stencil.m due to the legacy and convenience of coding 
% FC{16} The coordinates of stencil points of current face based on the face normal a 2_by_3 matrix = [coordinates of SP in downwind cell, coordinates of SP in upwind cell, coordinates of SP in further upwind cell]
% FC{17} The Zone ID of stencil points of current face based on the face normal
% a row vector = [Zone ID of SP in downwind cell, Zone ID of SP in upwind cell, Zone ID of SP in further upwind cell]
% FC{18} The three points that circle the stencil point a 4_by_3 matrix = [three points for SP in downwind cell, three points for SP in upwind cell, three points for SP in further upwind cell]
% For Each column vector,
% [ S1=The total number of boundary nodes in the three points:0-No node;1-one node;2-two nodes
%   S2=The order number of first point. If S1=0-order number of centroid; Else-order number of first boundary node
%   S3=The order number of second point. If S1<=1-order number of centroid; Else-order number of second boundary node
%   S4=The order number of third point. Always order of centroid, since S1<=2
% ]
% If S1=0, and S2=S3=S4, means the stencil is located at a centroid

%%%%%%%% FC{19}~{22} will be filled in stencil.m
% FC{19}   The assorted info of stencil points of current face based on V1
% FC{20}   The assorted info of stencil points of current face based on V2
%% Each of them is a S=cell(30,1) and each element in S is a row vectors follows:

%{S{1}  downwind cell # at each direction of V1 or V2; 
% S{2}  upwind cell # at each direction of V1 or V2; 
% S{3}  further upwind cell # at each direction of V1 or V2; 

% S{4}  The x coordinates of stencil point in downwind cell at each direction of V1 or V2;
% S{5}  The y coordinates of stencil point in downwind cell at each direction of V1 or V2;
% S{6}  The x coordinates of stencil point in upwind cell at each direction of V1 or V2;
% S{7}  The y coordinates of stencil point in upwind cell at each direction of V1 or V2;
% S{8}  The x coordinates of stencil point in further upwind cell at each direction of V1 or V2;
% S{9}  The y coordinates of stencil point in further upwind cell at each direction of V1 or V2;

% S{10} The zone ID of stencil point in downwind cell at each direction of V1 or V2;
% S{11} The zone ID of stencil point in upwind cell at each direction of V1 or V2;
% S{12} The zone ID of stencil point in further upwind cell at each direction of V1 or V2;

% S{13}  distance from the SP in downwind cell to the SP in upwind cell at each direction of V1 or V2;
% S{14}  distance from the SP in upwind cell to the SP in further upwind cell at each direction of V1 or V2; 
% S{15}  One over S{13};
% S{16}  One over S{14}; 

% S{17}  spatial ratio for current face,R1=(downwind cell size + upwind cell size)/upwind cell size,      at each direction of V1 or V2;
% S{18}  spatial ratio for extrapolation, R2=distance from the SP in downwind cell to the SP in upwind cell/ distance from the SP in upwind cell to the SP in further upwind cell, at each direction of V1 or V2;
% S{19}  One over S{17};
% S{20}  One over S{18};

% S{21} the left boundary node # on edge for extrapolation, at each direction of V1 or V2; 
% S{22} the right boundary node # on edge for extrapolation, at each direction of V1 or V2;
% S{23} spatial ratio for on edge,R3=distance from left node to intercept node/distance from left node to right node], at each direction of V1 or V2}

% FC{21}  The upwind identifier for 2nd-order mapping at the stencil points for V1
% a q1_by_3 matrix (q1 the number of velocities in V1) = [upwind identifier for SP in downwind cell, upwind identifier for SP in upwind cell, upwind identifier for SP in further upwind cell]
% For Each column vector, each entry is either 0 or 1
% 0---The SP is located upwind of the centroid, use 2nd-order mapping for SP
% 1---The SP is located downwind of the centroid, use the value at centroid for SP
% FC{22}  The upwind identifier for 2nd-order mapping at the stencil points for V2
% a q2_by_3 matrix (q2 the number of velocities in V2) = [upwind identifier for SP in downwind cell, upwind identifier for SP in upwind cell, upwind identifier for SP in further upwind cell]

%%%%%%%% FC{23}~{25} will be filled here
% FC{23}  The coordinates of stencil points of current face based on the face normal if one or many of the stencil points is out of boundaries
% FC{24}  The Zone ID of stencil points of current face based on the face normal if one or many of the stencil points is out of boundaries
% FC{25}  The three points that circle the stencil point if one or many of the stencil points is out of boundaries

%%%%%%%% FC{26}~{29} will be in stencil.m
% FC{26}  The assorted info of stencil points of current face based on V1 if one or many of the stencil points is out of boundaries
% FC{27}  The assorted info of stencil points of current face based on V2 if one or many of the stencil points is out of boundaries
% FC{28}  The upwind identifier for 2nd-order mapping at the stencil points for V1 if one or many of the stencil points is out of boundaries
% FC{29}  The upwind identifier for 2nd-order mapping at the stencil points for V2 if one or many of the stencil points is out of boundaries


%% Code starts
M=length(CELL);
N=length(NODE);
P=CELL{1};
e=ceil(10*sqrt(M*P{6})); %%%% Calculate the reference size
E=10^(numel(num2str((dx+dy)/20))-1);
O=length(edge(:,1));
%% Create the FACE data structure and fill in basic info from FC{1}~FC{11}
% Identify each distictive face, the idea is to use the order number of end
% nodes of each face, which is available in CELL
face_flag_2D=zeros(2,3*M); % Total number of possible faces, including duplicates
face_counter=0;
for i=1:M
    P=CELL{i};
    for k=1:3
        face_counter=face_counter+1;
        face_flag_2D(:,face_counter)=P{27+k};
    end
end
if face_counter~=3*M
    error('Wrong face counter!');
end

% Combine the two-dimension value into single-dimension in order to make
% the search easy, the idea is that since each face mid-point coordinate is unique, a
% number that combine combines the x and y coordinates will also be
% unique. It is not x+y, but x followed by y
face_flag_1D=zeros(1,3*M); % Total number of possible faces, including duplicates
for i=1:3*M
    face_flag_1D(i)=face_flag_2D(1,i)*E+face_flag_2D(2,i);
end
face_flag_1D_uniq=unique(face_flag_1D);
% Check uniqueness, this section is very slow, can be commented once the
% code is verified
if O~=length(face_flag_1D_uniq)
    error('Wrong number of faces!');
end
face_marker=zeros(O,1); % Contained integer values, 1 is boundary face, 2 is interior face, other values are error indicators
for i=1:M
    P=CELL{i};
    for k=1:3
        fc_id=P{27+k};
        face_id=fc_id(1)*E+fc_id(2);
        for j=1:O
            if face_id==face_flag_1D_uniq(j)
                break;
            end
        end
        if j==O
            if face_id~=face_flag_1D_uniq(j)
                error('face is not found!');
            end
        end
        face_marker(j)=face_marker(j)+1;
    end
end
face_marker_uniq=unique(face_marker);
if length(face_marker_uniq)==2
    if (face_marker_uniq(1)~=1 && face_marker_uniq(2)~=2) && (face_marker_uniq(1)~=2 && face_marker_uniq(2)~=1)
        error('The values should be 1 and 2');
    end
else
    error('There should be only two unique values!');
end

% Create FACE and fill FC{1}~FC{11}
K=55;
FACE=cell(K,1);
FACE={FACE,FACE};
for i=1:O
    FC=cell(K,1);
    FC{1,1}=i;
    face_found=0;
    for j=1:M
        P=CELL{j};
        for k=1:3
            fc_id=P{27+k};
            face_id=fc_id(1)*E+fc_id(2);
            if face_id==face_flag_1D_uniq(i)
                face_found=1;
                break;
            end
        end
        if face_found
            break;
        end
    end
    % Check exit condition
    if j==M && k==3
        if face_id~=face_flag_1D_uniq(i)
            error('Match face is not found!');
        end
    end
    % FC{5, 6, 8, 9, 10, 11}
    if k==1
        FC{5,1}=P{13};
        FC{6,1}=P{14};
        FC{8,1}=P{7};
        FC{9,1}=P{8};
        FC{10,1}=P{10};
        FC{11,1}=P{11};
    elseif k==2
        FC{5,1}=P{14};
        FC{6,1}=P{15};
        FC{8,1}=P{8};
        FC{9,1}=P{9};
        FC{10,1}=P{11};
        FC{11,1}=P{12};
    elseif k==3
        FC{5,1}=P{15};
        FC{6,1}=P{13};
        FC{8,1}=P{9};
        FC{9,1}=P{7};
        FC{10,1}=P{12};
        FC{11,1}=P{10};
    else
        error('There are only three faces!');
    end
    % FC{7}
    FC{7,1}=P{27+k};
    % FC{2}
    if FC{10,1}==1 && FC{11,1}==1
        mid_c=FC{7};
        if (mid_c(1)<X2 && mid_c(1)>X1) &&  (mid_c(2)<Y2 && mid_c(2)>Y1)
            FC{2,1}=0;
        else
            FC{2,1}=1;
        end
    elseif FC{10,1}==-1 && FC{11,1}==-1
        FC{2,1}=-1;
    elseif (FC{10,1}==-1 && FC{11,1}==1) || (FC{10,1}==1 && FC{11,1}==-1)
        FC{2,1}=0;
        disp('Check the meshing: the distance between the inner and outer boundaries may be too small!')
    else
        FC{2,1}=0;
    end
    % FC{3, 4}
    [Le,Ne]=en(FC{5,1},FC{6,1});
    FC{3,1}=Le;
    FC{4,1}=Ne;
    
    % Fill
    FACE{i}=FC;
end

%% Check FACE so far, using NODE
face_inner_counter=0;
face_interior_counter=0;
face_outer_counter=0;
for i=1:O
    FC=FACE{i};
    ND1=NODE{FC{8}};
    ND2=NODE{FC{9}};
    if ND1{2}~=FC{10} || ND2{2}~=FC{11}
        error('Region identifier of the end nodes of current face is wrong!');
    end
    if single(e+sum(ND1{3}==FC{5}))~=single(e+2) || single(e+sum(ND2{3}==FC{6}))~=single(e+2)
        error('Coordinates of the end nodes of current face are wrong!');
    end
    
    if single(e+sum(((ND1{3}+ND2{3})/2)==FC{7}))~=single(e+2)
        error('Coordinates of the mid-point of current face are wrong!');
    end
    
    if FC{2}==-1
        face_inner_counter=face_inner_counter+1;
    elseif FC{2}==0
        face_interior_counter=face_interior_counter+1;
    elseif FC{2}==1
        face_outer_counter=face_outer_counter+1;
    else
        error('false face region identifier!');
    end
end
if (face_inner_counter~=N_I_N || face_outer_counter~=(N-N_I)) || (face_inner_counter+face_interior_counter+face_outer_counter)~=O
    error('The face region identifier has false info!');
end

%% Fill FC{12}~FC{13}
for l=1:O
    FC=FACE{l};
    FC_NEIGH_CELL=zeros(2,2);
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
    if L==1
        Common_cell(end+1)=0; % On boundary
    end
    for i=1:L
        face_found_counter=0;
        P=CELL{Common_cell(i)};
        FC_MID1=P{28};
        FC_MID2=P{29};
        FC_MID3=P{30};
        %% Fill the first position in FC{12}, FC{13}
        Face_normal=-(P{5}'-FC{7}');
        if Face_normal*FC{4}'>0
            FC_NEIGH_CELL(1,1)=Common_cell(i);
            FC_NEIGH_CELL(2,1)=setxor(Common_cell(i),Common_cell);
        else
            FC_NEIGH_CELL(2,1)=Common_cell(i);
            FC_NEIGH_CELL(1,1)=setxor(Common_cell(i),Common_cell);
        end
        
        %% Fill the second position in FC{12}, FC{13}
        % Face 1
        if single(e+sum(FC_MID1==FC{7}))==single(e+2) % Changing
            face_found_counter=face_found_counter+1;
            if Face_normal*FC{4}'>0
                FC_NEIGH_CELL(1,2)=1;
            else
                FC_NEIGH_CELL(2,2)=1;
            end
        end
        % Face 2
        if single(e+sum(FC_MID2==FC{7}))==single(e+2) % Changing
            face_found_counter=face_found_counter+1;
            if Face_normal*FC{4}'>0
                FC_NEIGH_CELL(1,2)=2;
            else
                FC_NEIGH_CELL(2,2)=2;
            end
        end
        % Face 3
        if single(e+sum(FC_MID3==FC{7}))==single(e+2) % Changing
            face_found_counter=face_found_counter+1;
            if Face_normal*FC{4}'>0
                FC_NEIGH_CELL(1,2)=3;
            else
                FC_NEIGH_CELL(2,2)=3;
            end
        end
        if face_found_counter~=1
            error('The current face shared multiple times by the same cell, or the current face has no cell sharing it!');
        end
    end
    FC{12}=FC_NEIGH_CELL(1,:);
    FC{13}=FC_NEIGH_CELL(2,:);
    FACE{l}=FC;
end
% Check FC{12}~FC{13}
for l=1:O
    FC=FACE{l};
    neigh_up=FC{12};
    neigh_down=FC{13};
    if length(union(neigh_up(1,1),neigh_down(1,1)))==2
        if length(union(0,[neigh_up(1,1),neigh_down(1,1)]))==2 % boundary face, inner boundary or outer boundary
            cell_id=setxor(0,[neigh_up(1,1),neigh_down(1,1)]);
            face_ord=setxor(0,[neigh_up(1,2),neigh_down(1,2)]);
            P=CELL{cell_id};
            if P{1+face_ord}~=0
                error('This face should be on boundary!');
            end
        elseif length(union(0,[neigh_up(1,1),neigh_down(1,1)]))==3 % interior face
            P_up=CELL{neigh_up(1,1)};
            P_down=CELL{neigh_down(1,1)};
            if P_up{1+neigh_up(1,2)}~=P_down{1} || P_down{1+neigh_down(1,2)}~=P_up{1}
                error('Wrong pair!');
            end
            n_c=P_down{5}-P_up{5};
            if FC{4}*n_c<=0
                error('The face normal should be reversed!');
            end
        else
            error('Logic error!');
        end
    else
        error('The neighor cells of the same face should have different IDs!');
    end
end


%% Fill FC{14}~FC{15}
for l=1:O
    FC=FACE{l};
    if FC{2}==1
        if FC{10}~=1 || FC{11}~=1
            error('The boundary identifier of the face is not consistent with its two end nodes!');
        end
        FC_NEIGH_CELL1=FC{12};
        FC_NEIGH_CELL2=FC{13};
        ND1=NODE{FC{8}};
        ND2=NODE{FC{9}};
        cell_star_1=ND1{9};
        cell_star_2=ND2{9};
        common_cell=intersect(cell_star_1,cell_star_2);
        if length(common_cell)~=2
            error('The face on periodic boundary should have two attached cells!');
        end
        if FC_NEIGH_CELL1(1,1)==0 && FC_NEIGH_CELL2(1,1)~=0
            if length(setxor(common_cell,FC_NEIGH_CELL2(1,1)))~=1
                error('Wrong neighbor is found!');
            end
            FC{14}=setxor(common_cell,FC_NEIGH_CELL2(1,1));
            P=CELL{FC{14}};
            for i=1:3
                DH=abs(P{27+i}-FC{7});
                if (single(e+DH(1))==single(e) && single(e+DH(2))==single(e+(Y2-Y1))) || (single(e+DH(1))==single(e+(X2-X1)) && single(e+DH(2))==single(e))
                    break;
                end
            end
            if i==3
                if (single(e+DH(1))~=single(e) && single(e+DH(2))~=single(e+(Y2-Y1))) && (single(e+DH(1))~=single(e+(X2-X1)) && single(e+DH(2))~=single(e))
                    error('The face is not on boundary!');
                end
            end
            FC{14}=[FC{14},i];
            
            FC{15}=FC{13};
        elseif FC_NEIGH_CELL1(1,1)~=0 && FC_NEIGH_CELL2(1,1)==0
            if length(setxor(common_cell,FC_NEIGH_CELL1(1,1)))~=1
                error('Wrong neighbor is found!');
            end
            FC{14}=FC{12};
            FC{15}=setxor(common_cell,FC_NEIGH_CELL1(1,1));
            P=CELL{FC{15}};
            for i=1:3
                DH=abs(P{27+i}-FC{7});
                if (single(e+DH(1))==single(e) && single(e+DH(2))==single(e+(Y2-Y1))) || (single(e+DH(1))==single(e+(X2-X1)) && single(e+DH(2))==single(e))
                    break;
                end
            end
            if i==3
                if (single(e+DH(1))~=single(e) && single(e+DH(2))~=single(e+(Y2-Y1))) && (single(e+DH(1))~=single(e+(X2-X1)) && single(e+DH(2))~=single(e))
                    error('The face is not on boundary!');
                end
            end
            FC{15}=[FC{15},i];
        else
            error('The input for FC{12} and FC{13} or current boundary face is incorrect!');
        end
        FACE{l}=FC;
    end
end
% Check FC{14}~FC{15}
for l=1:O
    FC=FACE{l};
    if FC{2}==1 
        neigh_up=FC{14};
        neigh_down=FC{15};
        P_up=CELL{neigh_up(1,1)};
        P_down=CELL{neigh_down(1,1)};
        DH=abs(P_up{27+neigh_up(1,2)}-P_down{27+neigh_down(1,2)});
        if (single(e+DH(1))~=single(e) && single(e+DH(2))~=single(e+(Y2-Y1))) && (single(e+DH(1))~=single(e+(X2-X1)) && single(e+DH(2))~=single(e))
            error('The periodic face neighbor is not correct!');
        end
    end
end