function [ele,neigh,Star]=PostTri_4Tri_killer(ele,node,neigh,Star)
% function [ele,neigh,Star]=PostTri_4Tri_killer(ele,node,neigh,Star) kills
% all the structure that a node with four connected triangles by swapping
% one of the 4-triangle with its neighbor
% ele--the ele file from Triangle
% node---the node file from Triangle
% neigh--the neigh file from triangle
% Star---the cell structure data for all star cell structure

M=length(ele(:,1));
N=length(node(:,1));
FourTri_counter=0;
for r=1:N
    P=Star{r};
    if P{1,1}==4 && node(r,4)==0
        FourTri_counter=FourTri_counter+1;
        nd_4Tri(FourTri_counter)=r;
    end
end
%% Plot mesh & highlight the nodes that have 4 triangles connected
edge1_x=zeros(M,2);
edge2_x=zeros(M,2);
edge3_x=zeros(M,2);
edge1_y=zeros(M,2);
edge2_y=zeros(M,2);
edge3_y=zeros(M,2);
for r=1:M
    nd1=ele(r,2);
    nd2=ele(r,3);
    nd3=ele(r,4);
    edge1_x(r,:)=[node(nd1,2),node(nd2,2)];
    edge2_x(r,:)=[node(nd2,2),node(nd3,2)];
    edge3_x(r,:)=[node(nd3,2),node(nd1,2)];
    edge1_y(r,:)=[node(nd1,3),node(nd2,3)];
    edge2_y(r,:)=[node(nd2,3),node(nd3,3)];
    edge3_y(r,:)=[node(nd3,3),node(nd1,3)];
end
figure;
for r=1:M
    plot(edge1_x(r,:),edge1_y(r,:));
    hold on
    plot(edge2_x(r,:),edge2_y(r,:));
    hold on
    plot(edge3_x(r,:),edge3_y(r,:));
    hold on
end
hold on
for r=1:FourTri_counter
    plot(node(nd_4Tri(r),2),node(nd_4Tri(r),3), 'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'red','MarkerEdgeColor', 'red');
    hold on
end
axis equal tight
% Check whether any two or three found nodes belong to the same triangle
if FourTri_counter>=2
    for r=1:FourTri_counter
        P=Star{nd_4Tri(r)};
        star=P{2,1};
        for i=1:P{1,1}
            T=ele(star(i),:); % Each of the triangles connected to the node that has 4 connected triangles
            dup_nd_4Tri=0;
            for d=1:3
                for s=1:FourTri_counter
                    if T(1,2+d-1)==nd_4Tri(s)
                        dup_nd_4Tri=dup_nd_4Tri+1;
                    end
                end
            end
            if dup_nd_4Tri>=2
                msg=['The node ', num2str(nd_4Tri(r)), ' and at least another node has seperately 4 connected triangles and these nodes belong to the same triangle.']
                error('Remesh with Triangle!');
            end
        end
    end
end
% Check whether any of total 4 triangles connected to a node is the neighbor of a triangle that is any of another
% total 4 triangles connected to another node
pair_counter=0;
if FourTri_counter>=2
    for r=1:FourTri_counter
        P=Star{nd_4Tri(r)};
        star=P{2,1};
        for i=1:P{1,1}
            T=ele(star(i),:); % Each of the triangles connected to the node that has 4 connected triangles
            count=0;
            pair_found=0;
            for d=1:3 % Find the esge with two nodes whose neighbor is going to be found
                if T(1,2+d-1)~=nd_4Tri(r)
                    count=count+1;
                    edge_lookup_neigh(count)=T(1,2+d-1);
                end
            end
            if length(edge_lookup_neigh)~=2
                error('The edge for looking for neighbor is not found!');
            end
            % with all the found nodes that have 4 connected triangles, to
            % find which node has one of its triangle the neighbor of
            % current edge
            for s=1:FourTri_counter
                if s~=r % Exclude the current node that has 4 connected triangles
                    Q=Star{nd_4Tri(s)};
                    star_q=Q{2,1};
                    for j=1:Q{1,1}
                        T_q=ele(star_q(j),:);
                        count_q=0;
                        for d=1:3 % Find the esge with two nodes whose neighbor is going to be found
                            if T_q(1,2+d-1)~=nd_4Tri(s)
                                count_q=count_q+1;
                                edge_found_neigh(count_q)=T_q(1,2+d-1);
                            end
                        end
                        if length(edge_found_neigh)~=2
                            error('The edge of the found neighbor is not found!');
                        end
                        if (edge_lookup_neigh(1)==edge_found_neigh(1) && edge_lookup_neigh(2)==edge_found_neigh(2)) || (edge_lookup_neigh(2)==edge_found_neigh(1) && edge_lookup_neigh(1)==edge_found_neigh(2))
                            pair_counter=pair_counter+1;
                            nd_4Tri_del_pair(pair_counter,1)=nd_4Tri(r);
                            nd_4Tri_del_pair(pair_counter,2)=nd_4Tri(s);
                            Tri_swap_pair(pair_counter,1)=star(i);
                            Tri_swap_pair(pair_counter,2)=star_q(j);
                            pair_found=1;
                            break;
                        end
                    end
                    if pair_found
                        break;
                    end
                end
            end
        end
    end
end
% Delete the duplicate pair
if pair_counter>1
    c_t=pair_counter;
    for i=1:c_t
        nd_4Tri_del_pair_temp=nd_4Tri_del_pair(i+1:end,:);
        NN=length(nd_4Tri_del_pair_temp(:,1));
        c=0;
        bi=[Inf Inf];
        for j=1:NN
            if norm(nd_4Tri_del_pair(i,:))~=norm(nd_4Tri_del_pair_temp(j,:))
                c=c+1;
                bi(c,:)=nd_4Tri_del_pair_temp(j,:);
            end
        end
        if sum(bi==Inf)==2
            if i==1
                nd_4Tri_del_pair=nd_4Tri_del_pair(i,:);
            else
                ;
            end
        else
            nd_4Tri_del_pair=[nd_4Tri_del_pair(1:i,:);bi];
        end
        c_t=length(nd_4Tri_del_pair(:,1));
    end
    
    c_t=pair_counter;
    for i=1:c_t
        Tri_swap_pair_temp=Tri_swap_pair(i+1:end,:);
        NN=length(Tri_swap_pair_temp(:,1));
        c=0;
        bi=[Inf Inf];
        for j=1:NN
            if norm(Tri_swap_pair(i,:))~=norm(Tri_swap_pair_temp(j,:))
                c=c+1;
                bi(c,:)=Tri_swap_pair_temp(j,:);
            end
        end
        if sum(bi==Inf)==2
            if i==1
                Tri_swap_pair=Tri_swap_pair(i,:);
            else
                ;
            end
        else
            Tri_swap_pair=[Tri_swap_pair(1:i,:);bi];
        end
        c_t=length(Tri_swap_pair(:,1));
    end
    % check
    pair_counter_1=length(nd_4Tri_del_pair(:,1));
    pair_counter_2=length(Tri_swap_pair(:,1));
    if pair_counter_1~=pair_counter_2
        error('Pairs are not matched!');
    else
        pair_counter=pair_counter_1;
    end
end
% Do swap and update nd_4Tri & FourTri_counter
if pair_counter>0
    % Swap all the triangles pairs
    for i=1:pair_counter
        ne1=Tri_swap_pair(i,1);
        ne2=Tri_swap_pair(i,2);
        [ele,neigh]=PostTri_swap(ele,neigh,ne1,ne2);
    end
    % change the nd_4Tri_del_pair into a row vector and sort out all unique
    % node
    nd_4Tri_del=unique(nd_4Tri_del_pair);
    % delete those unique nodes from nd_4Tri
    nd_4Tri_del_counter=length(nd_4Tri_del);
    loc_del=zeros(1,FourTri_counter);
    for r=1:FourTri_counter
        for t=1:nd_4Tri_del_counter
            if nd_4Tri(r)==nd_4Tri_del(t)
                loc_del(r)=r;
                break;
            end
        end
    end
    FourTri_counter_new=0;
    for r=1:FourTri_counter
        if loc_del(r)~=r
            FourTri_counter_new=FourTri_counter_new+1;
            nd_4Tri_new(FourTri_counter_new)=nd_4Tri(r);
        end
    end
    FourTri_counter=FourTri_counter_new;
    nd_4Tri=nd_4Tri_new;
end
% Check if there exists a triangle who is the neighbor of two triangles that belong to two seperate 4Tri structure
% No swap is performed, just mark up which triangle pair should not be
% swapped
not_swap=zeros(M,1);
for i=1:FourTri_counter
    Nd=nd_4Tri(i);
    P=Star{Nd};
    if P{1,1}~=4
        error('The current node does not have 4 connected triangles!');
    end
    star=P{2,1};
    for j=1:4
        cn=0;
        T=ele(star(j),:);
        which_edge=zeros(1,3);
        for k=1:3
            if T(1+k)~=Nd
                cn=cn+1;
                edge_nd(cn)=T(1+k);
            else
                which_edge(k)=1; % which edge is going to be swapped
            end
        end
        if cn~=2
            error('The edge for swapping is not found!');
        end
        for k=1:3
            if which_edge(k)==1
                e_num=k;
                break;
            end
        end
        neighbor=neigh(star(j),e_num+1);
        % Check whether this neighbor has a neighbor that is one of the
        % triangle of another 4-Tri structure
        commTri_neigh_counter=0;
        for k=1:3
            if neigh(neighbor,1+k)~=star(j)
                commTri_neigh_counter=commTri_neigh_counter+1;
                commTri_neigh(commTri_neigh_counter)=neigh(neighbor,1+k);
            end
        end
        if commTri_neigh_counter~=2
            error('The neighbor of common triangle should exclude the current triangle!');
        end
        for k=1:commTri_neigh_counter
            for ii=1:FourTri_counter
                if i~=ii % Exclude the current 4_tri structure
                    ND=nd_4Tri(ii);
                    Q=Star{ND};
                    star_Q=Q{2,1};
                    for jj=1:Q{1,1}
                        if star_Q(jj)==commTri_neigh(k) % The neighbor of common triangle is member of another 4-tri structure 
                            not_swap(star(j))=1;
                            break;
                        end
                    end
                end
                if not_swap(star(j))
                    break;
                end
            end
            if not_swap(star(j))
                break;
            end
        end
    end
end
% Every node that has 4 connected triangles are far away enough from each other then the regular swap scheme can be applied
% including the scenario that there exists a triangle who is the neighbor of two triangles that belong to two seperate 4Tri structure
% No swap is performed, just mark up which triangle pair should not be swapped
for i=1:FourTri_counter
    Nd=nd_4Tri(i);
    P=Star{Nd};
    if P{1,1}~=4
        error('The current node does not have 4 connected triangles!');
    end
    star=P{2,1};
    for j=1:4
        cn=0;
        T=ele(star(j),:);
        which_edge=zeros(1,3);
        for k=1:3
            if T(1+k)~=Nd
                cn=cn+1;
                edge_nd(cn)=T(1+k);
            else
                which_edge(k)=1; % which edge is going to be swapped
            end
        end
        if cn~=2
            error('The edge for swapping is not found!');
        end
        nd1=node(edge_nd(1),:);
        nd2=node(edge_nd(2),:);
        Q1=Star{edge_nd(1)};
        Q2=Star{edge_nd(2)};
        S1=Q1{4};
        S2=Q2{4};
        sat_nd_on_boundary=0;
        for s=1:Q1{3}
            if node(S1(s),4)~=0
                sat_nd_on_boundary=1;
                break;
            end
        end
        for s=1:Q2{3}
            if node(S2(s),4)~=0
                sat_nd_on_boundary=1;
                break;
            end
        end
        if (Q1{1}>5 && Q2{1}>5) && (nd1(4)==0 || nd2(4)==0) && (~not_swap(star(j))) && (~sat_nd_on_boundary)
            for k=1:3
                if which_edge(k)==1
                    e_num=k;
                    break;
                end
            end
            neighbor=neigh(star(j),e_num+1); % first colum in neigh file is ele order number
            % check if the found neighbor triangle is one of the 4-triangle
            % structure
            for s=1:4
                if neighbor==star(s)
                    error('The neighbor can not be one of the 4-tri strucuture!');
                end
            end
            [ele,neigh]=PostTri_swap(ele,neigh,star(j),neighbor);
            break;
        end
    end
end
% Update Star
for r=1:N
    P=Star{r};
    [n,star]=PostTri_findstar(r,ele);
    P{1,1}=n; % number of triangles in the current star structure
    P{2,1}=star; % vector of those triangle order numbers
    q=0;
    poly_star=0;
    for l=1:n
        T=ele(star(l),2:4);
        for w=1:3
            if T(w)~=r
                if l==1 && (w==1 || w==2)
                    q=q+1;
                    poly_star(q)=T(w);
                else
                    found_dup_nd=0;
                    for x=1:q
                        if T(w)==poly_star(x)
                            found_dup_nd=1;
                        end
                    end
                    if found_dup_nd
                        ;
                    else
                        q=q+1;
                        poly_star(q)=T(w);
                    end
                end
            end
        end
    end
    P{3,1}=q; % The number of saterlite nodes in the current star structure
    P{4,1}=poly_star; % The vector of saterlite node order numbers in the current star structure
    
    Star{r}=P;
end
% Check
FourTri_counter=0;
for r=1:N
    P=Star{r};
    if P{1,1}==4 && node(r,4)==0
        FourTri_counter=FourTri_counter+1;
    end
end
if FourTri_counter~=0
    error('Still exists 4-triangle structure');
end