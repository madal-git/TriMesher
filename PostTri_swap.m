function [ele,neigh]=PostTri_swap(ele,neigh,ne1,ne2)
% function [ele,neigh]=PostTri_swap(ele,neigh,ne1,ne2) swap the element ne1 and ne2
% in the ele file from Triangle and and update ele and neigh file.
% ele is the ele file from Triangle
% neigh is the neigh file from Triangle
% ne1 is an interger, the order number of first triangle to be swapped
% ne2 is an interger, the order number of second triangle to be swapped

% Check if two different triangles are given to swap
if ne1==ne2
    error('Can not swap single triangle!');
end
% Check whether ne1 and ne2 are neighbors and on which edge they are
% neighbors
T1=ele(ne1,:);
T2=ele(ne2,:);
N1=neigh(ne1,:);
N2=neigh(ne2,:);
neighbor_found=0;
for i=1:3
    % The edge of first triangle
    if i==1
        edge1=[T1(1,3),T1(1,4)];
    elseif i==2
        edge1=[T1(1,4),T1(1,2)];
    else
        edge1=[T1(1,2),T1(1,3)];
    end
    for j=1:3
        % The edge of second triangle
        if j==1
            edge2=[T2(1,3),T2(1,4)];
        elseif j==2
            edge2=[T2(1,4),T2(1,2)];
        else
            edge2=[T2(1,2),T2(1,3)];
        end
        % Check if two edges are the same edge, if yes, two triangles are
        % neighbors
        if edge1(1)==edge2(2) && edge1(2)==edge2(1) % It will never happen that edge1(1)==edge2(1) && edge1(2)==edge2(2)
            % because the nodes in both triangles are numbered
            % counterclockwisely
            neighbor_found=1;
            break;
        else
            ;
        end
    end
    if neighbor_found
        break;
    end
end
if neighbor_found
    % Node i of T1 and Node j of T2 will form the new shared edge of new two neughbor triangles
    new_edge_nd1=T1(1,i+1);
    new_edge_nd2=T2(1,j+1);
    T1_new=T1;
    T2_new=T2;
    N1_new=N1;
    N2_new=N2;
    % Create two new triangles, both of their first node is new_edge_nd1
    if i==1
        % First new triangle
        T1_new(1,2)=new_edge_nd1;
        T1_new(1,3)=T1(1,3);
        T1_new(1,4)=new_edge_nd2;
        % Second new triangle
        T2_new(1,2)=new_edge_nd1;
        T2_new(1,3)=new_edge_nd2;
        T2_new(1,4)=T1(1,4);
        % neighbor
        if j==1
            N1_new(1,2)=N2(1,3);
            N1_new(1,3)=ne2;
            N1_new(1,4)=N1(1,4);
            
            N2_new(1,2)=N2(1,4);
            N2_new(1,3)=N1(1,3);
            N2_new(1,4)=ne1;
        elseif j==2
            N1_new(1,2)=N2(1,4);
            N1_new(1,3)=ne2;
            N1_new(1,4)=N1(1,4);
            
            N2_new(1,2)=N2(1,2);
            N2_new(1,3)=N1(1,3);
            N2_new(1,4)=ne1;
        else
            N1_new(1,2)=N2(1,2);
            N1_new(1,3)=ne2;
            N1_new(1,4)=N1(1,4);
            
            N2_new(1,2)=N2(1,3);
            N2_new(1,3)=N1(1,3);
            N2_new(1,4)=ne1;
        end
    elseif i==2
        % First new triangle
        T1_new(1,2)=new_edge_nd1;
        T1_new(1,3)=T1(1,4);
        T1_new(1,4)=new_edge_nd2;
        % Second new triangle
        T2_new(1,2)=new_edge_nd1;
        T2_new(1,3)=new_edge_nd2;
        T2_new(1,4)=T1(1,2);
        % neighbor
        if j==1
            N1_new(1,2)=N2(1,3);
            N1_new(1,3)=ne2;
            N1_new(1,4)=N1(1,2);
            
            N2_new(1,2)=N2(1,4);
            N2_new(1,3)=N1(1,4);
            N2_new(1,4)=ne1;
        elseif j==2
            N1_new(1,2)=N2(1,4);
            N1_new(1,3)=ne2;
            N1_new(1,4)=N1(1,2);
            
            N2_new(1,2)=N2(1,2);
            N2_new(1,3)=N1(1,4);
            N2_new(1,4)=ne1;
        else
            N1_new(1,2)=N2(1,2);
            N1_new(1,3)=ne2;
            N1_new(1,4)=N1(1,2);
            
            N2_new(1,2)=N2(1,3);
            N2_new(1,3)=N1(1,4);
            N2_new(1,4)=ne1;
        end
    else
        % First new triangle
        T1_new(1,2)=new_edge_nd1;
        T1_new(1,3)=T1(1,2);
        T1_new(1,4)=new_edge_nd2;
        % Second new triangle
        T2_new(1,2)=new_edge_nd1;
        T2_new(1,3)=new_edge_nd2;
        T2_new(1,4)=T1(1,3);
        % neighbor
        if j==1
            N1_new(1,2)=N2(1,3);
            N1_new(1,3)=ne2;
            N1_new(1,4)=N1(1,3);
            
            N2_new(1,2)=N2(1,4);
            N2_new(1,3)=N1(1,2);
            N2_new(1,4)=ne1;
        elseif j==2
            N1_new(1,2)=N2(1,4);
            N1_new(1,3)=ne2;
            N1_new(1,4)=N1(1,3);
            
            N2_new(1,2)=N2(1,2);
            N2_new(1,3)=N1(1,2);
            N2_new(1,4)=ne1;
        else
            N1_new(1,2)=N2(1,2);
            N1_new(1,3)=ne2;
            N1_new(1,4)=N1(1,3);
            
            N2_new(1,2)=N2(1,3);
            N2_new(1,3)=N1(1,2);
            N2_new(1,4)=ne1;
        end
    end
else
   error('The given two triangles are not neighbors. Can not swap!');
end
%% The other two neighbor update of T1
if i==1
    N2neigh_T1_1=neigh(N1(1,4),:);
    N2neigh_T1_2=neigh(N1(1,3),:);
elseif i==2
    N2neigh_T1_1=neigh(N1(1,2),:);
    N2neigh_T1_2=neigh(N1(1,4),:);
else
    N2neigh_T1_1=neigh(N1(1,3),:);
    N2neigh_T1_2=neigh(N1(1,2),:);
end
N2neigh_T1_found_1=0;
for ii=1:3
    if N2neigh_T1_1(1,ii+1)==ne1
        N2neigh_T1_1(1,ii+1)=ne1;
        N2neigh_T1_found_1=1;
        break;
    end
end
N2neigh_T1_found_2=0;
for ii=1:3
    if N2neigh_T1_2(1,ii+1)==ne1
        N2neigh_T1_2(1,ii+1)=ne2;
        N2neigh_T1_found_2=1;
        break;
    end
end
if N2neigh_T1_found_1==0 || N2neigh_T1_found_2==0
    error('The other two neighbors of T1 is not found!');
end
neigh(N2neigh_T1_1(1,1),:)=N2neigh_T1_1;
neigh(N2neigh_T1_2(1,1),:)=N2neigh_T1_2;
%% The other two neighbor update of T2
if j==1
    N2neigh_T2_1=neigh(N2(1,3),:);
    N2neigh_T2_2=neigh(N2(1,4),:);
elseif j==2
    N2neigh_T2_1=neigh(N2(1,4),:);
    N2neigh_T2_2=neigh(N2(1,2),:);
else
    N2neigh_T2_1=neigh(N2(1,2),:);
    N2neigh_T2_2=neigh(N2(1,3),:);
end
N2neigh_T2_found_1=0;
for jj=1:3
    if N2neigh_T2_1(1,jj+1)==ne2
        N2neigh_T2_1(1,jj+1)=ne1;
        N2neigh_T2_found_1=1;
        break;
    end
end
N2neigh_T2_found_2=0;
for jj=1:3
    if N2neigh_T2_2(1,jj+1)==ne2
        N2neigh_T2_2(1,jj+1)=ne2;
        N2neigh_T2_found_2=1;
        break;
    end
end
if N2neigh_T2_found_1==0 || N2neigh_T2_found_2==0
    error('The other two neighbors of T1 is not found!');
end
neigh(N2neigh_T2_1(1,1),:)=N2neigh_T2_1;
neigh(N2neigh_T2_2(1,1),:)=N2neigh_T2_2;
%% Update ele file & neigh file
ele(ne1,:)=T1_new;
ele(ne2,:)=T2_new;
neigh(ne1,:)=N1_new;
neigh(ne2,:)=N2_new;