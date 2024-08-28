function [element,M,node,NIB,NI,N,neighbor]=PostTri_move_rec(element,node,neighbor)
% function [element,M,node,NIB,NI,N,neighbor]=PostTri_move_rec(element,node,neighbor)
% moves the interior nodes in the data sequence before outer boundary nodes
% but after inner boundary nodes if there are holes. Then change the value
% and data sequence in element file and neighbor file accordingly. 
% In addition, M, NIB, NI and N is generated as the variables for the solver.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% node %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Integrety check of node file
[N,L]=size(node);
% Total number of inner boundary nodes
inner=0;
for r=1:N
    nd=node(r,:);
    if nd(1,L)==-1
        inner=inner+1;
    end
end
% Total number of interior boundary nodes
interior=0;
for r=1:N
    nd=node(r,:);
    if nd(1,L)==0
        interior=interior+1;
    end
end
% Total number of outer boundary nodes
outer=0;
for r=1:N
    nd=node(r,:);
    if nd(1,L)==1
        outer=outer+1;
    end
end
% Check total number
if inner+interior+outer~=N
    error('Check the boundary marker of node!');
end


%% Sample nodes for each category
% Inner boundary nodes
node_inn=zeros(inner,L);
a=0;
for r=1:N
    nd=node(r,:);
    if nd(1,L)==-1
        a=a+1;
        node_inn(a,:)=nd;
    end
end
if a~=inner
    error('sampling is not complete!');
end
% Interior nodes
node_int=zeros(interior,L);
b=0;
for r=1:N
    nd=node(r,:);
    if nd(1,L)==0
        b=b+1;
        node_int(b,:)=nd;
    end
end
if b~=interior
    error('sampling is not complete!');
end
% Outer boundary nodes
node_out=zeros(outer,L);
c=0;
for r=1:N
    nd=node(r,:);
    if nd(1,L)==1
        c=c+1;
        node_out(c,:)=nd;
    end
end
if c~=outer
    error('sampling is not complete!');
end
%% Move and assembly
if inner==0 % no holes
    node_new=[node_int;node_out];
else
    node_new=[node_inn;node_int;node_out];
end
% Creat the mapping of each node's old position and new position
map_nd=zeros(N,2);
for r=1:N
    nd=node_new(r,:);
    map_nd(r,1)=nd(1,1); % Old position
    map_nd(r,2)=r; % New position
end
%% Creat the new node file
for r=1:N
    node_new(r,1)=r;
end
node=node_new;
NIB=inner;
NI=inner+interior;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% node %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% element %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% chang the node number in the element file
[M,S]=size(element);
for r=1:M
    elt=element(r,:);
    for k=1:S
        if k==1
            ;
        else
            for t=1:N
                if elt(1,k)==map_nd(t,1)
                    elt(1,k)=map_nd(t,2);
                    break;
                end
            end
            % Check whether matching is found or just the loop ends
            if t==N
                if elt(1,k)~=map_nd(t,2)
                    error('The matching node is not found in the current element!');
                end
            end
        end
    end
    element(r,:)=elt;
end
%% make the arrangment of three vertices of each element from counterclockwise to clockwise
% Here choose to switch first and second vortex
temp=element(:,2);
element(:,2)=element(:,3);
element(:,3)=temp;
%% make arrangment accordingly to neighbor file
temp=neighbor(:,2);
neighbor(:,2)=neighbor(:,3);
neighbor(:,3)=temp;