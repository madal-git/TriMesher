function [buff_layer_nd,buff_layer_seg,N]=PreTri_buff_inner(base_nd,FBuf_inner,HL,Origin)
% function [buff_layer_nd]=Pretri_buff_inner(base_nd) create a layer of interior
% nodes outside of the base_nd. buff_layer_nd and base_nd has the same data
% structure, where each row stands for each node, [col1, col2, col3, col4].
% col1----The node order number
% col2----x
% col3----y
% col4----boundary marker. -1----Inner boundary; 0----Interior node
% base_nd is the given base nodes outside of which the buffer layler will
% be created
% FBuf_inner is the total number of buffer layers in the inner hole
% HL is the paramters about the inner hole
% Origin is the origin of the mesh domain
% buff_layer_nd is the created node file of the buffer layer
% buff_layer_seg is the created segment file of the buffer layer


%% Check legitimacy of base_nd
[N,L]=size(base_nd);
if L~=4
    error('The data structure of given base nodes is not recoganizable!');
end
%% Generate the buffer layer
% Initialize
buff_layer_nd=zeros(N,L);
buff_layer_seg=zeros(N,L);
buff_vector=zeros(N,L); % The vector parameters for each node, [col1, col2, col3, col4].
% col1---The x component of normal direction toward which the a buffer
% layer node should be located
% col2---The y component of normal direction toward which the a buffer
% layer node should be located
% col3---The x component of tangential direction toward which the a buffer
% layer node should be moved
% col4---The y component of tangential direction toward which the a buffer
% layer node should be moved
buff_dis=zeros(N,2); % The distance that determines the location of the buffer layer nodes
% [col1; col2] col1---The distance along the normal direction starting from
% the base node; col2---The distance the the buffer layer node should be
% moved along the tangential direction

% Determine the parameters buff_vector and buff_dis
for i=1:N
    if i==1
        nd_up=base_nd(N,2:3);
        nd_down=base_nd(i+1,2:3);
    elseif i==N
        nd_up=base_nd(i-1,2:3);
        nd_down=base_nd(1,2:3);
    else
        nd_up=base_nd(i-1,2:3);
        nd_down=base_nd(i+1,2:3);
    end
    nd_center=base_nd(i,2:3);
    
    nn_x_updown=-(nd_up(2)-nd_down(2))/dis(nd_up',nd_down');
    nn_y_updown=(nd_up(1)-nd_down(1))/dis(nd_up',nd_down');
    nt_x_updown=(nd_down(1)-nd_up(1))/dis(nd_up',nd_down');
    nt_y_updown=(nd_down(2)-nd_up(2))/dis(nd_up',nd_down');
    
    nn_x_up=-(nd_up(2)-nd_center(2))/dis(nd_up',nd_center');
    nn_y_up=(nd_up(1)-nd_center(1))/dis(nd_up',nd_center');
    nt_x_up=(nd_center(1)-nd_up(1))/dis(nd_up',nd_center');
    nt_y_up=(nd_center(2)-nd_up(2))/dis(nd_up',nd_center');
    
    nn_x_down=-(nd_center(2)-nd_down(2))/dis(nd_center',nd_down');
    nn_y_down=(nd_center(1)-nd_down(1))/dis(nd_center',nd_down');
    nt_x_down=(nd_down(1)-nd_center(1))/dis(nd_center',nd_down');
    nt_y_down=(nd_down(2)-nd_center(2))/dis(nd_center',nd_down');
    
%     nn_x=(nn_x_updown+nn_x_up+nn_x_down)/3;
%     nn_y=(nn_y_updown+nn_y_up+nn_y_down)/3;
%     nt_x=(nt_x_updown+nt_x_up+nt_x_down)/3;
%     nt_y=(nt_y_updown+nt_y_up+nt_y_down)/3;
    % Alternative proach
    nn_x=(nn_x_up+nn_x_down)/2;
    nn_y=(nn_y_up+nn_y_down)/2;
    nt_x=(nt_x_up+nt_x_down)/2;
    nt_y=(nt_y_up+nt_y_down)/2;
    
    buff_vector(i,:)=[nn_x,nn_y,nt_x,nt_y];
    dis_seg=dis(nd_down',[base_nd(i,2);base_nd(i,3)]);
    dis_n=dis_seg*sqrt(3)/2;
    dis_t=dis_seg/2;
    buff_dis(i,:)=[dis_n,dis_t];
end
% Create buffer layer node one by one. Two steps, 1. create it along the
% normal direction; 2. move along tangential direction.
if FBuf_inner<40
    for i=1:N
        buff_layer_nd(i,1)=i;
        temp_nd(:,1)=[base_nd(i,2);base_nd(i,3)]+buff_dis(i,1)*[buff_vector(i,1);buff_vector(i,2)]; % Step 1
        temp_nd(:,1)=temp_nd(:,1)+buff_dis(i,2)*[buff_vector(i,3);buff_vector(i,4)]; % Step 2
        buff_layer_nd(i,2)=temp_nd(1,1);
        buff_layer_nd(i,3)=temp_nd(2,1);
        buff_layer_nd(i,4)=0; % Marker for interior node
    end
else % When the # of buffer layers is more than 40, the numerical errors of 
     % calculating moving direction and distance will diverge, therefore alternative
     % is required
     if HL{1}==0 % Rectangular hole
         ;
     elseif HL{1}==1 % Circular hole
         % Check the whether the hole violates the outer boundaries
         Ori_cir_x=Origin(1)+HL{3};
         Ori_cir_y=Origin(2)+HL{4};
         dis_seg=dis([base_nd(1,2);base_nd(1,3)],[base_nd(2,2);base_nd(2,3)]);
         dis_n=dis_seg*sqrt(3)/2;
         R=dis([Ori_cir_x;Ori_cir_y],[base_nd(1,2);base_nd(1,3)])+dis_n;
         Agl_ini=angle([(base_nd(1,2)-Ori_cir_x);(base_nd(1,3)-Ori_cir_y)],[1;0])/180*pi;
         dA=2*pi/N;
         Agl_ini=Agl_ini+dA/2;
         for i=1:N
             buff_layer_nd(i,1)=i;
             buff_layer_nd(i,2)=Ori_cir_x+R*cos(Agl_ini+(i-1)*dA);
             buff_layer_nd(i,3)=Ori_cir_y+R*sin(Agl_ini+(i-1)*dA);
             buff_layer_nd(i,4)=0; % Marker for interior node
         end
     else
         error('Other holes than rectangle or circle is currently not available!');
     end
end
% Create the segment file
for r=1:N
    buff_layer_seg(r,1)=r;
    if r==N
        buff_layer_seg(r,2)=buff_layer_nd(r,1);
        buff_layer_seg(r,3)=buff_layer_nd(1,1);
    else
        buff_layer_seg(r,2)=buff_layer_nd(r,1);
        buff_layer_seg(r,3)=buff_layer_nd(r+1,1);
    end
    buff_layer_seg(r,4)=0;
end