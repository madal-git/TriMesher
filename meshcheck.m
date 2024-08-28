function meshcheck(CELL,NODE,FACE,N_H,N_L,X,Y,FM)
% function meshcheck(CELL,NODE,FACE,N_H,RatioLH,X,Y,FM) check the legitimacy
% of generated mesh data, before delivering it to the solver.

% CELL is the cell data structure
% NODE is the node data structure
% FACE is the face data structure
% N_H is the total number of nodes on the vertical outer boundary
% N_H is the total number of nodes on the vertical outer boundary
% X and Y is vector that stores the x and y coordinates of all nodes.
% FM is the flag for mesh type. FM=0----IRT mesh; FM=1----General mesh

X1=min(X);
X2=max(X);
Y1=min(Y);
Y2=max(Y);
M=length(CELL);
N=length(NODE);
O=length(FACE);
P=CELL{1};
e=ceil(10*sqrt(M*P{6})); %%%% Calculate the reference size
e1=10;

if FM==0
    %% Checking integraty of the mesh
    Area_real=(X2-X1)*(Y2-Y1);
    Area_actual=0;
    for r=1:M
        P=CELL{r};
        Area_actual=Area_actual+P{6};
    end
    if single(e+Area_actual)~=single(e+Area_real)
        error('The integrity of the mesh domain is voilated!');
    end
    %% Checking the number of key elements
    % number of cells
    N_cell=(N_H-1)*(N_L-1)*4;
    if N_cell~=M
        error('The total number of cells is incorrect!');
    end
    % number of nodes
    N_node=N_H*N_L+(N_H-1)*(N_L-1);
    if N_node~=N
        error('The total number of nodes is incorrect!');
    end
    % number of faces
    N_face=(N_H-1)*(N_L-1)*4+(N_L-1)*N_H+(N_H-1)*N_L;
    if N_face~=O
        error('The total number of faces is incorrect!');
    end
    %% Check boundary identifiers
    % Cell data
    N_boundary_node=0;
    N_interior_node=0;
    N_boundary_face=0;
    N_interior_face=0;
    
    ID_node=0;
    ID_face=0;
    
    for r=1:M
        P=CELL{r};
        for i=1:3
            node_b_id=P{9+i};
            face_b_id=P{18+i};
            if node_b_id==0
                ID_node(end+1)=P{6+i};
                if length(ID_node)==length(unique(ID_node))
                    N_interior_node=N_interior_node+1;
                end
                ID_node=unique(ID_node);
            elseif node_b_id==1;
                ID_node(end+1)=P{6+i};
                if length(ID_node)==length(unique(ID_node))
                    N_boundary_node=N_boundary_node+1;
                end
                ID_node=unique(ID_node);
            else
                error('The boundary identifier could ONLY be 0 or 1!');
            end
            % 
            if face_b_id==0
                ID_face(end+1)=P{15+i};
                if length(ID_face)==length(unique(ID_face))
                    N_interior_face=N_interior_face+1;
                end
                ID_face=unique(ID_face);
            elseif face_b_id==1;
                ID_face(end+1)=P{15+i};
                if length(ID_face)==length(unique(ID_face))
                    N_boundary_face=N_boundary_face+1;
                end
                ID_face=unique(ID_face);
            else
                error('The boundary identifier could ONLY be 0 or 1!');
            end
        end
    end
    if N_boundary_node~=(2*(N_H+N_L)-4) || (N_boundary_node+N_interior_node)~=N
        error('The boundary identifier for nodes has false information!');
    end
    if N_boundary_face~=(2*(N_H+N_L-2)) || (N_boundary_face+N_interior_face)~=O
        error('The boundary identifier for facs has false information!');
    end
    % Cell data
    N_boundary_node=0;
    N_interior_node=0;
    for s=1:N
        ND=NODE{s};
        if ND{2}==0
            N_interior_node=N_interior_node+1;
        elseif ND{2}==1;
            N_boundary_node=N_boundary_node+1;
        else
            error('The boundary identifier could ONLY be 0 or 1!');
        end
    end
    if N_boundary_node~=(2*(N_H+N_L)-4) || (N_boundary_node+N_interior_node)~=N
        error('The boundary identifier for nodes has false information!');
    end
    % Face data
    N_boundary_node=0;
    N_interior_node=0;
    N_boundary_face=0;
    N_interior_face=0;
    
    ID_node=0;
    
    for t=1:O
        FC=FACE{t};
        if FC{2}==0
            N_interior_face=N_interior_face+1;
        elseif FC{2}==1
            N_boundary_face=N_boundary_face+1;
        else
            error('The boundary identifier could ONLY be 0 or 1!');
        end
        for i=1:2
            node_b_id=FC{9+i};
            if node_b_id==0
                ID_node(end+1)=FC{7+i};
                if length(ID_node)==length(unique(ID_node))
                    N_interior_node=N_interior_node+1;
                end
                ID_node=unique(ID_node);
            elseif node_b_id==1;
                ID_node(end+1)=FC{7+i};
                if length(ID_node)==length(unique(ID_node))
                    N_boundary_node=N_boundary_node+1;
                end
                ID_node=unique(ID_node);
            else
                error('The boundary identifier could ONLY be 0 or 1!');
            end
        end
    end
    if N_boundary_node~=(2*(N_H+N_L)-4) || (N_boundary_node+N_interior_node)~=N
        error('The boundary identifier for nodes has false information!');
    end
    if N_boundary_face~=(2*(N_H+N_L-2)) || (N_boundary_face+N_interior_face)~=O
        error('The boundary identifier for facs has false information!');
    end
    %% Check connections of nodal facial information throughout all data structures
    for r=1:M
        P=CELL{r};
        % Nodal info
        for i=1:3
            fc_id=P{15+i};
            FC=FACE{fc_id};
            nd1_id=FC{8};
            nd2_id=FC{9};
            ND1=NODE{nd1_id};
            ND2=NODE{nd2_id};
            if ND1{2}~=FC{10} || ND2{2}~=FC{11}
                error('The boundary identifier for the same node in NODE and FACE do not match!');
            end
            if i==1
                if P{7}==ND1{1} && P{8}==ND2{1}
                    if P{10}~=ND1{2} || P{11}~=ND2{2}
                        error('The boundary identifier for the same node in NODE and CELL do not match!');
                    end
                    Bool_coord1=single(e+P{13})~=single(e+ND1{3});
                    Bool_coord2=single(e+P{14})~=single(e+ND2{3});
                    if (Bool_coord1(1) || Bool_coord1(end)) || (Bool_coord2(1) || Bool_coord2(end))
                        error('The coordinates for the same face in CELL and FACE do not match!');
                    end
                elseif P{8}==ND1{1} && P{7}==ND2{1}
                    if P{11}~=ND1{2} || P{12}~=ND2{2}
                        error('The boundary identifier for the same node in NODE and CELL do not match!');
                    end
                    Bool_coord1=single(e+P{14})~=single(e+ND1{3});
                    Bool_coord2=single(e+P{13})~=single(e+ND2{3});
                    if (Bool_coord1(1) || Bool_coord1(end)) || (Bool_coord2(1) || Bool_coord2(end))
                        error('The coordinates for the same face in CELL and FACE do not match!');
                    end
                else
                    error('The order number for the same node in NODE and CELL do not match!');
                end
            elseif i==2
                if P{8}==ND1{1} && P{9}==ND2{1}
                    if P{11}~=ND1{2} || P{12}~=ND2{2}
                        error('The boundary identifier for the same node in NODE and CELL do not match!');
                    end
                    Bool_coord1=single(e+P{14})~=single(e+ND1{3});
                    Bool_coord2=single(e+P{15})~=single(e+ND2{3});
                    if (Bool_coord1(1) || Bool_coord1(end)) || (Bool_coord2(1) || Bool_coord2(end))
                        error('The coordinates for the same face in CELL and FACE do not match!');
                    end
                elseif P{9}==ND1{1} && P{8}==ND2{1} % the order of two ends in FC is not corresponding to the order of node in P
                    if P{12}~=ND1{2} || P{11}~=ND2{2}
                        error('The boundary identifier for the same node in NODE and CELL do not match!');
                    end
                    Bool_coord1=single(e+P{15})~=single(e+ND1{3});
                    Bool_coord2=single(e+P{14})~=single(e+ND2{3});
                    if (Bool_coord1(1) || Bool_coord1(end)) || (Bool_coord2(1) || Bool_coord2(end))
                        error('The coordinates for the same face in CELL and FACE do not match!');
                    end
                else
                    error('The order number for the same node in NODE and CELL do not match!');
                end
            else
                if P{9}==ND1{1} && P{7}==ND2{1}
                    if P{12}~=ND1{2} || P{10}~=ND2{2}
                        error('The boundary identifier for the same node in NODE and CELL do not match!');
                    end
                    Bool_coord1=single(e+P{15})~=single(e+ND1{3});
                    Bool_coord2=single(e+P{13})~=single(e+ND2{3});
                    if (Bool_coord1(1) || Bool_coord1(end)) || (Bool_coord2(1) || Bool_coord2(end))
                        error('The coordinates for the same face in CELL and FACE do not match!');
                    end
                elseif P{7}==ND1{1} && P{9}==ND2{1} % the order of two ends in FC is not corresponding to the order of node in P
                    if P{10}~=ND1{2} || P{12}~=ND2{2}
                        error('The boundary identifier for the same node in NODE and CELL do not match!');
                    end
                    Bool_coord1=single(e+P{13})~=single(e+ND1{3});
                    Bool_coord2=single(e+P{15})~=single(e+ND2{3});
                    if (Bool_coord1(1) || Bool_coord1(end)) || (Bool_coord2(1) || Bool_coord2(end))
                        error('The coordinates for the same face in CELL and FACE do not match!');
                    end
                else
                    error('The order number for the same node in NODE and CELL do not match!');
                end
            end
        end
        % Facial info
        for i=1:3
            fc_id=P{15+i};
            FC=FACE{fc_id};
            if FC{2}~=P{18+i}
                error('The boundary identifier for the same face in CELL and FACE do not match!');
            end
            Bool_normal=single(e+abs(FC{4}))~=single(e+abs(P{33+i}));
            Bool_mid=single(e+FC{7})~=single(e+P{27+i});
            if single(e+FC{3})~=single(e+P{30+i}) || (Bool_normal(1) || Bool_normal(end)) || (Bool_mid(1) || Bool_mid(end))
                error('The geometric info for the same face in CELL and FACE do not match!');
            end
        end
    end
    %% Check centroid, area, face length and normal for CV
    for r=1:M
        P=CELL{r};
        bool_cent=(single(e+P{5})~=single(e+(P{13}+P{14}+P{15})/3));
        if bool_cent(1) || bool_cent(end)
            error('The centroid of current cell is not correct!');
        end
        if single(e+P{6})~=single(e+area_tri(P{13},P{14},P{15}))
            error('The area of current cell is not correct!');
        end
        for i=1:3
            if i==1
                coor_nd1=P{13};
                coor_nd2=P{14};
            elseif i==2
                coor_nd1=P{14};
                coor_nd2=P{15};
            else
                coor_nd1=P{15};
                coor_nd2=P{13};
            end
            if single(e+dis(coor_nd1,coor_nd2))~=single(e+P{30+i})
                error('The lenght of current face is incorrect!');
            end
            n_out=P{27+i}'-P{5}';
            if n_out*P{33+i}'<0 || single(e+P{33+i}*(coor_nd1-coor_nd2))~=single(e)
                error('The direction of current face normal is incorrect!');
            end
            if single(e+norm(P{33+i},2))~=single(e+1)
                error('The magnitude of current face normal is incorrect!');
            end
        end
    end
    %% Check neighbor information in the 1st column of FC{12}~FC{15}
    for r=1:M
        P=CELL{r};
        for i=1:3
            FC=FACE{P{15+i}};
            if P{1+i}==0
                neighbor1=FC{14};
                neighbor2=FC{15};
                neigh_pair_face1=[neighbor1(1,1),neighbor2(1,1)];
                if length(unique(neigh_pair_face1))~=2 || length(unique(union(neigh_pair_face1,0)))==2
                    error('error in two neighbors of current face!');
                end
                cell_p_neigh_id=setxor(neigh_pair_face1,P{1});
                if length(cell_p_neigh_id)~=1
                    error('Periodic neighbor is not found!');
                end
                C=CELL{cell_p_neigh_id};
                for j=1:3
                    if C{1+j}==0
                        break;
                    end
                end
                if j==3
                    if C{1+j}~=0
                        error('Periodic neighbor is not found!');
                    end
                end
                FC_p=FACE{C{15+j}};
                neighbor_p1=FC_p{14};
                neighbor_p2=FC_p{15};
                neigh_pair_face2=[neighbor_p1(1,1),neighbor_p2(1,1)];
                if length(unique(neigh_pair_face2))~=2 || length(unique(union(neigh_pair_face2,0)))==2
                    error('error in two neighbors of current face!');
                end
                if length(unique(union(neigh_pair_face1,neigh_pair_face2)))~=2
                    error('The face on the boundary and its corresponding face do not have the same pair of neighbor cells!');
                end
            else
                neigh_pair_cell=[P{1},P{1+i}];
                neighbor1=FC{12};
                neighbor2=FC{13};
                neigh_pair_face=[neighbor1(1,1),neighbor2(1,1)];
                if length(unique(neigh_pair_cell))~=length(unique(union(neigh_pair_cell,neigh_pair_face)))
                    error('The neighbor info in CELL and FACE does not match!');
                end
            end
        end
    end
    %% Check neighbor information in the 2nd column of FC{12}~FC{15}, and the sequencing of 1st & 2nd neighbor cell of current face
    for r=1:O
        FC=FACE{r};
        if FC{2}==0
            neighbor1=FC{12};
            neighbor2=FC{13};
            cell1=CELL{neighbor1(1,1)};
            cell2=CELL{neighbor2(1,1)};
            face1=FACE{cell1{15+neighbor1(1,2)}};
            face2=FACE{cell2{15+neighbor2(1,2)}};
            if face1{1}~=FC{1} || face2{1}~=FC{1}
                error('The info about ith face of the pair neighbor cell is incorrect!');
            end
            face1_normal=cell1{33+neighbor1(1,2)};
            face2_normal=cell2{33+neighbor2(1,2)};
            bool=single(e+face1_normal+face2_normal)~=single(e);
            if bool(1) || bool(end)
                error('The info about ith face of the pair neighbor cell is incorrect!');
            end
            if double(e1+face1_normal*FC{4}')~=double(e1+1) || double(e1+face2_normal*FC{4}')~=double(e1-1)
                error('The 1st and 2nd neighbor cells of current face should be switched!');
            end
        elseif FC{2}==1
            neighbor1=FC{14};
            neighbor2=FC{15};
            cell1=CELL{neighbor1(1,1)};
            cell2=CELL{neighbor2(1,1)};
            face1=FACE{cell1{15+neighbor1(1,2)}};
            face2=FACE{cell2{15+neighbor2(1,2)}};
            mid1=face1{7};
            mid2=face2{7};
            if single(e+dis(mid1,mid2))~=single(e+(X2-X1)) && single(e+dis(mid1,mid2))~=single(e+(Y2-Y1))
                error('The info about ith face of the pair neighbor cell on the boundary is incorrect!');
            end
            face1_normal=cell1{33+neighbor1(1,2)};
            face2_normal=cell2{33+neighbor2(1,2)};
            bool=single(e+face1_normal+face2_normal)~=single(e);
            if bool(1) || bool(end)
                error('The info about ith face of the pair neighbor cell is incorrect!');
            end
            if single(e1+face1_normal*FC{4}')~=single(e1+1) || single(e1+face2_normal*FC{4}')~=single(e1-1)
                error('The 1st and 2nd neighbor cells of current face should be switched!');
            end
        else
            error('The boundary ID for face could only be 0 or 1 for the time being!');
        end
    end
    %% Check star structure, do not need to check the periodic boundary node
    for s=1:N
        ND=NODE{s};
        Star_face=ND{13};
        Star_cell=0;
        for i=1:ND{12}
            FC=FACE{Star_face(i)};
            neighbor1=FC{12};
            neighbor2=FC{13};
            Star_cell=union(Star_cell,[neighbor1(1,1),neighbor2(1,1)]);
        end
        Star_cell=setxor(unique(Star_cell),0);
        if length(unique(ND{5}))~=length(unique(union(ND{5},Star_cell)))
            error('The star structure for the current interior node has false info!');
        end
    end
elseif FM==1
    ;
else
    error('The flag for mesh type is incorrect!');
end