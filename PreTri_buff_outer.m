function [buff_layer_nd,buff_layer_seg,K]=PreTri_buff_outer(base_nd)
% [buff_layer_nd,buff_layer_seg,N]=PreTri_buff_outer(base_nd) create a layer of interior
% nodes inside of the base_nd that are on the rectangular outer boundary. 
% buff_layer_nd and base_nd has the same data
% structure, where each row stands for each node, [col1, col2, col3, col4].
% col1----The node order number
% col2----x
% col3----y
% col4----boundary marker. -1----Inner boundary; 0----Interior node
% base_nd is the given base nodes inside of which the buffer layler will
% be created
% buff_layer_nd is the created node file of the buffer layer
% buff_layer_seg is the created segment file of the buffer layer


%% Check legitimacy of data structure of base_nd
[S,D]=size(base_nd);
if D~=4
    error('The data structure of given base nodes is not recoganizable!');
end
%% Check whether base_nd forms a rectangular closed loop
for r=1:S
    X(r)=base_nd(r,2);
    Y(r)=base_nd(r,3);
end
X1=min(X);
X2=max(X);
Y1=min(Y);
Y2=max(Y);
c_top=0;
c_right=0;
c_bottom=0;
c_left=0;
for r=1:S
    if single(base_nd(r,3))==single(Y2)
        c_top=c_top+1;
    end
    if single(base_nd(r,2))==single(X2)
        c_right=c_right+1;
    end
    if single(base_nd(r,3))==single(Y1)
        c_bottom=c_bottom+1;
    end
    if single(base_nd(r,2))==single(X1)
        c_left=c_left+1;
    end
end
if (c_top+c_right+c_bottom+c_left-4)~=S
    error('The base nodes do not form a rectangle1');
end
% Check spacing
if c_top~=c_bottom
    error('The top and bottom surface have different spacing!');
else
    NL=c_top;
end

if c_left~=c_right
    error('The left and right surface have different spacing!');
else
    NH=c_left;
end
Dx=(X2-X1)/(NL-1);
Dy=(Y2-Y1)/(NH-1);
%% Generate the buffer layer Node file
X1_new=X1+Dy*sqrt(3)/2;
X2_new=X2-Dy*sqrt(3)/2;
Y1_new=Y1+Dx*sqrt(3)/2;
Y2_new=Y2-Dx*sqrt(3)/2;
NL_new=NL-1;
NH_new=NH-1;
Dx_new=(X2_new-X1_new)/(NL_new-1);
Dy_new=(Y2_new-Y1_new)/(NH_new-1);
K=(NL_new+NH_new-2)*2;
if K~=(S-4)
    error('The total number of nodes on the buffer layer is incorrect!');
end
buff_layer_nd=zeros(K,D);
Top_start=(base_nd(1,2)+base_nd(2,2))/2;
Right_start=(base_nd(NL,3)+base_nd(NL+1,3))/2;
Bottom_start=(base_nd(NL+NH-1,2)+base_nd(NL+NH,2))/2;
Left_start=(base_nd(NL+NH-1+NL-1,3)+base_nd(NL+NH-1+NL,3))/2;
for r=1:K
    buff_layer_nd(r,1)=r;
    if r<=NL_new-2
        buff_layer_nd(r,2)=X1_new+r*Dx_new;
        buff_layer_nd(r,3)=Y2_new;
    elseif r==NL_new-1
        buff_layer_nd(r,2)=X2_new;
        buff_layer_nd(r,3)=Y2_new;
    elseif r<=NL_new-1+NH_new-2
        buff_layer_nd(r,2)=X2_new;
        buff_layer_nd(r,3)=Y2_new-(r-NL_new+1)*Dy_new;
    elseif r==NL_new-1+NH_new-1
        buff_layer_nd(r,2)=X2_new;
        buff_layer_nd(r,3)=Y1_new;
    elseif r<=NL_new-1+NH_new-1+NL_new-2
        buff_layer_nd(r,2)=X2_new-(r-NL_new+1-NH_new+1)*Dx_new;
        buff_layer_nd(r,3)=Y1_new;
    elseif r==NL_new-1+NH_new-1+NL_new-1
        buff_layer_nd(r,2)=X1_new;
        buff_layer_nd(r,3)=Y1_new;
    elseif r<K
        buff_layer_nd(r,2)=X1_new;
        buff_layer_nd(r,3)=Y1_new+(r-NL_new+1-NH_new+1-NL_new+1)*Dy_new;
    else
        buff_layer_nd(r,2)=X1_new;
        buff_layer_nd(r,3)=Y2_new;
    end
    buff_layer_nd(r,4)=0;
%     figure(1);
%     plot(buff_layer_nd(r,2),buff_layer_nd(r,3),'Marker', 'o','Markersize',4, 'MarkerFaceColor', 'green','MarkerEdgeColor', 'green');
%     axis equal tight;
%     hold on;
end
% Create the buffer layer segment file
buff_layer_seg=zeros(K,D);
for r=1:K
    buff_layer_seg(r,1)=r;
    if r==K
        buff_layer_seg(r,2)=buff_layer_nd(r,1);
        buff_layer_seg(r,3)=buff_layer_nd(1,1);
    else
        buff_layer_seg(r,2)=buff_layer_nd(r,1);
        buff_layer_seg(r,3)=buff_layer_nd(r+1,1);
    end
    buff_layer_seg(r,4)=0;
end