function F=PostTri_cell_rec(element,node,neighbor,M)
% function F=PostTri_cell_rec(element,node,neighbor) create the cell data F
% for FVLBM solver
% element---the processed element data
% node---the processed node data
% neighbor---the processed neighbor data
% M---the total number of triangles

K=40;
F=cell(K,1);
F={F,F};

for r=1:M
    P=cell(K,1);
    P{1,1}=element(r,1);
    nd1=element(r,2);
    nd2=element(r,3);
    nd3=element(r,4);
    
    P{7,1}=nd1;
    P{10,1}=node(nd1,4);
    P{13,1}=[node(nd1,2);node(nd1,3)];
    
    P{8,1}=nd2;
    P{11,1}=node(nd2,4);
    P{14,1}=[node(nd2,2);node(nd2,3)];
    
    P{9,1}=nd3;
    P{12,1}=node(nd3,4);
    P{15,1}=[node(nd3,2);node(nd3,3)];
    
    x1=[node(nd1,2),node(nd2,2)];
    x2=[node(nd2,2),node(nd3,2)];
    x3=[node(nd3,2),node(nd1,2)];
    y1=[node(nd1,3),node(nd2,3)];
    y2=[node(nd2,3),node(nd3,3)];
    y3=[node(nd3,3),node(nd1,3)];
    
    P{22,1}=x1;   % x coordinate of face 1
    P{23,1}=y1;   % y coordinate of face 1
    
    P{24,1}=x2;   % x coordinate of face 2
    P{25,1}=y2;   % y coordinate of face 2
    
    P{26,1}=x3;   % x coordinate of face 3
    P{27,1}=y3;   % y coordinate of face 3
    
    P{5,1}=(P{13,1}+P{14,1}+P{15,1})/3; % Centroid coordinate of Cell 1
    
    P{28,1}=(P{13,1}+P{14,1})/2;    % Mid point coordinate of face 1
    P{29,1}=(P{14,1}+P{15,1})/2;    % Mid point coordinate of face 2
    P{30,1}=(P{15,1}+P{13,1})/2;    % Mid point coordinate of face 3
    
    P{6,1}=tri_area(P{13,1},P{14,1},P{15,1});   % area of cell
    
    % neighbor number on each face
    % The order starts as 4-2-3, instead of 2-3-4 is because that in
    % Triangle, the first neighbor cell is defined as the cell opposing the
    % first vertex. However, in this code, the first neighbor cell is
    % defined as the cell attached to the first face that is constructed by
    % connecting the first and second vertices. So the order 4-2-3 will
    % rotate the defination.
    P{2,1}=neighbor(r,4);
    if P{2,1}==-1
        P{2,1}=0;
    end
    P{3,1}=neighbor(r,2);
    if P{3,1}==-1
        P{3,1}=0;
    end
    P{4,1}=neighbor(r,3);
    if P{4,1}==-1
        P{4,1}=0;
    end
    
    F{r}=P;
end