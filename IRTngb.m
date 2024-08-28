function [CELL,FACE]=IRTngb(CELL,NODE,FACE)
% [CELL,FACE]=IRTngb(CELL,NODE,FACE) find the neighbor information and fill it
% into the cell and face data structures.
% CELL is the cell structure data of all triangles.
% NODE is the node structure data of all triangles.
% FACE is the face structure data of all triangles.

%%% The data to be filled is listed as follows.
%%%%%%%%%% In CELL, C=CELL{m}
% C{2}    The order number of neighbor triangle on face 1(0:No neighbor; nonezero:has neighbor)
% C{3}    The order number of neighbor triangle on face 2(0:No neighbor; nonezero:has neighbor)
% C{4}    The order number of neighbor triangle on face 3(0:No neighbor; nonezero:has neighbor)

%%%%%%%%%% In FACE, FC=FACE{o}
% FC{12}  The order # of 1st neighbor cell (0:No neighbor; nonezero:has neighbor)
% FC{13}  The order # of 2nd neighbor cell (0:No neighbor; nonezero:has neighbor)
%%%%%%%%%%%%%%%%%%%%%%%%%% For Periodic boundaries ONLY %%%%%%%%%%%%%%%%%%%%%%%%%%%
% FC{14}  The order # of 1st neighbor cell (empty for interior face;nonezero on boundary and FC{2} has to be 1)
% FC{15}  The order # of 2nd neighbor cell (empty for interior face;nonezero on boundary and FC{2} has to be 1)

M=length(CELL);
N=length(NODE);
O=length(FACE);
P=CELL{1};
e=ceil(10*sqrt(M*P{6})); %%%% Calculate the reference size

%% Fill C{2}~C{4}, and FC{12}~FC{13}
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
        FC_ID1=P{16};
        FC_ID2=P{17};
        FC_ID3=P{18};
        %% Fill the first position in FC{12}, FC{13}
        Face_normal=-(P{5}'-FC{7}');
        if Face_normal*FC{4}'>0
            FC_NEIGH_CELL(1,1)=Common_cell(i);
            FC_NEIGH_CELL(2,1)=setxor(Common_cell(i),Common_cell);
        else
            FC_NEIGH_CELL(2,1)=Common_cell(i);
            FC_NEIGH_CELL(1,1)=setxor(Common_cell(i),Common_cell);
        end
        
        %% Fill the second position in FC{12}, FC{13}, and C{2}, C{3}, C{4}
        % Face 1
        if FC{1}==FC_ID1 % Changing
            face_found_counter=face_found_counter+1;
            if Face_normal*FC{4}'>0
                FC_NEIGH_CELL(1,2)=1;
            else
                FC_NEIGH_CELL(2,2)=1;
            end
            P{2}=setxor(Common_cell(i),Common_cell); % Changing
        end
        % Face 2
        if FC{1}==FC_ID2 % Changing
            face_found_counter=face_found_counter+1;
            if Face_normal*FC{4}'>0
                FC_NEIGH_CELL(1,2)=2;
            else
                FC_NEIGH_CELL(2,2)=2;
            end
            P{3}=setxor(Common_cell(i),Common_cell); % Changing
        end
        % Face 3
        if FC{1}==FC_ID3 % Changing
            face_found_counter=face_found_counter+1;
            if Face_normal*FC{4}'>0
                FC_NEIGH_CELL(1,2)=3;
            else
                FC_NEIGH_CELL(2,2)=3;
            end
            P{4}=setxor(Common_cell(i),Common_cell); % Changing
        end
        if face_found_counter~=1
            error('The current face shared multiple times by the same cell, or the current face has no cell sharing it!');
        end
        CELL{Common_cell(i)}=P;
    end
    FC{12}=FC_NEIGH_CELL(1,:);
    FC{13}=FC_NEIGH_CELL(2,:);
    FACE{l}=FC;
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
                Face_normal=P{33+i};
                if single(e+abs(Face_normal*FC{4}'))==single(e) || single(e+abs(Face_normal*FC{4}'))==single(e+1)
                    break;
                end
            end
            if i==3
                Face_normal=P{33+i};
                if single(e+abs(Face_normal*FC{4}'))~=single(e) && single(e+abs(Face_normal*FC{4}'))~=single(e+1)
                    error('Periodic face is not found!');
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
                Face_normal=P{33+i};
                if single(e+abs(Face_normal*FC{4}'))==single(e) || single(e+abs(Face_normal*FC{4}'))==single(e+1)
                    break;
                end
            end
            if i==3
                Face_normal=P{33+i};
                if single(e+abs(Face_normal*FC{4}'))~=single(e) && single(e+abs(Face_normal*FC{4}'))~=single(e+1)
                    error('Periodic face is not found!');
                end
            end
            FC{15}=[FC{15},i];
        else
            error('The input for FC{12} and FC{13} or current boundary face is incorrect!');
        end
        FACE{l}=FC;
    end
end