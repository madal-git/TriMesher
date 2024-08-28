function [ele,neigh,Star]=PostTri_7Tri_killer(ele,node,neigh,Star)
% function [ele,neigh,Star]=PostTri_7Tri_killer(ele,node,neigh,Star) kills
% all the structure that a node with seven connected triangles by swapping
% one of the 7-triangle with its neighbor
% ele--the ele file from Triangle
% node---the node file from Triangle
% neigh--the neigh file from triangle
% Star---the cell structure data for all star cell structure

M=length(ele(:,1));
N=length(node(:,1));
SevenTri_counter=0;
nd_7Tri=0;
for r=1:N
    P=Star{r};
    if P{1,1}==7 && node(r,4)==0
        SevenTri_counter=SevenTri_counter+1;
        nd_7Tri(SevenTri_counter)=r;
    end
end
if length(nd_7Tri)>1
    SevenTri_found=1;
else
    if nd_7Tri~=0
        SevenTri_found=1;
    else
        SevenTri_found=0;
    end
end
SevenTri_counter
a=0;
% while SevenTri_found
    a=a+1
    %% Study all satalite nodes to see how many triangles are connected to each satalite node and find the satalite node Q that forms the swap pair with central node P
    P=Star{nd_7Tri(1)};
    poly_star=P{4,1};
%     SixTri_poly_counter=0;
%     nd_6Tri_poly=0;
    for i=1:P{3,1}
        Q=Star{poly_star(i)};
        if Q{1,1}==7 && node(poly_star(i),4)==0
            break;
        else
            if Q{1,1}==6 && node(poly_star(i),4)==0
%                 SixTri_poly_counter=SixTri_poly_counter+1;
%                 nd_6Tri_poly(SixTri_poly_counter)=poly_star(i);
               break;
            end
        end
    end
    common_nodes=[nd_7Tri(1),poly_star(i)];
    %% Find the two elements in the star strucure that shares points P and Q
    star=P{2,1};
    pair=0;
    pair_counter=0;
    for j=1:P{1,1}
        vertices=ele(star(j),2:4);
        if length(union(vertices,common_nodes))==3
            pair_counter=pair_counter+1;
            pair(pair_counter)=star(j);
        end
    end
    if length(pair)~=2
        error('Swapping pair is not found!');
    end
    ne1=pair(1);
    ne2=pair(2);
    [ele,neigh]=PostTri_swap(ele,neigh,ne1,ne2);
    %% Update Star
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
    %% Find 7-tri strucure again
    nd_7Tri=0;
    for r=1:N
        P=Star{r};
        if P{1,1}==7 && node(r,4)==0
            nd_7Tri=r;
            break;
        end
    end
    if nd_7Tri~=0
        SevenTri_found=1;
    else
        SevenTri_found=0;
    end
% end

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

    plot(node(nd_7Tri,2),node(nd_7Tri,3), 'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'red','MarkerEdgeColor', 'red');
    hold on
axis equal tight