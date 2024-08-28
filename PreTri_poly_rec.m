function [Nd, Seg, Hole]=PreTri_poly_rec(H,RatioL2H,N_H,N_L,Origin,Fd,Num_holes,holes_para,FBuf_inner,FBuf_outer,Fcob,Fw)

% function [Nd, Seg, Hole]=PreTri_poly_rec(H,RatioL2H,N_H,N_L,Origin,Fd,Num_holes,holes_para,FBuf_inner,FBuf_outer,Fcob,Fw)
% prepares the boundary node definition for the tool Triangle to generate
% the triangular unstructured mesh. 
% Nd is the matrix for all the boundary nodes information. Each row of the
% matrix is for each node. [node #, x, y, boundary marker]
% Seg is another matrix for all segments information. Each segment connect
% two nearest boundary nodes. Each row is for each segment. [seg #, nd1, nd2, boundary marker]
% Hole is matrix for information of all holes. Each row is for each hole.
% [hole #, x, y] where x and y are the coordinates of any point inside of
% the hole.


% H is the height of the rectangular channel.
% RatioL2H is the ratio of the length of channel over the height.
% N_H is the wished total number of nodes on the left or right surface of
% the channel.
% N_L is the wished total number of nodes on the top or bottom surface of
% the channel.
% Origin is the coordinate of the bottom-left corner of the channel.
% Num_holes is the total number of holes inside of the channel. If it is
% zero, the last variable will be ignored.
% Fd is the tag for tuning the spacing on the outer rectangular boundaries.
% holes_para is a cell data type that stores all information about each
% hole. It's a cell structure with each cell for each hole. [cell 1, cell2, cell 3, ...]. cell 1
% indicate the hole structure, 0---rectangular, 1---circular, others is tbd. The data
% structure from cell 2 for different hole structure will be different. For
% rectangular, [cell1, cell2---height, cell3---the ratio of length over height, cell4---the distance from
% the left surface of the hole to the left boundary of the channel, cell5---the distance from
% the bottom surface of the hole to the bottom boundary of the channel, cell6---wished total number of 
% nodes on left or right surface of the hole, cell7---wished total number of 
% nodes on top or bottom surface of the hole]
% For circular holes, [cell1, cell2---radius of the hole, cell3---the distance from
% the center of the hole to the left boundary of the channel, cell4---the distance from
% the center of the hole to the bottom boundary of the channel, cell5-wished number of nodes on the
% circular surface]


%% Check conditions of holes parameters
if Num_holes==0
    ;
else
    S=size(holes_para);
    if Num_holes~=S(1,2)
        error('The number of holes is more or less than what is allowed!');
    end
end
X1=Origin(1);
X2=X1+H*RatioL2H;
Y1=Origin(2);
Y2=Y1+H;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generating nodes and segments on holes
if Num_holes==0
    HoleNd=0;
    HoleSeg=0;
    Hole=0;
else
    % Initialize the counter
    Counter=0;
    for k=1:Num_holes
        if Num_holes==1
            HL=holes_para;
        else
            HL=holes_para(:,k);
        end
        
        if HL{1}==0 % Rectangular hole
            ;
        elseif HL{1}==1 % Circular hole
            % Check the whether the hole violates the outer boundaries
            Ori_cir_x=Origin(1)+HL{3};
            Ori_cir_y=Origin(2)+HL{4};
            if Ori_cir_x+HL{2}>=X2 || Ori_cir_x-HL{2}<=X1 || Ori_cir_y+HL{2}>=Y2 || Ori_cir_y-HL{2}<=Y1
                error('The current hole is beyond the outer boundaries! Move the hole or decrease the size of it!');
            end
            % Generate nodes counterclockwise on the circular hole
            %HoleNd=zeros(HL{5},4);
            dA=2*pi/HL{5}; % Angle spacing on the circular hole.
            Agl_ini=0; % The position where the first node is located.
            for r=1:HL{5}
                HoleNd(Counter+r,1)=r+Counter;
                HoleNd(Counter+r,2)=Ori_cir_x+HL{2}*cos(Agl_ini+(r-1)*dA);
                HoleNd(Counter+r,3)=Ori_cir_y+HL{2}*sin(Agl_ini+(r-1)*dA);
                HoleNd(Counter+r,4)=-1;
%                 figure(1);
%                 plot(HoleNd(Counter+r,2),HoleNd(Counter+r,3),'Marker', 'o','Markersize',4, 'MarkerFaceColor', 'blue');
%                 axis equal tight;
%                 hold on;
            end
            % Generate the segments
            %HoleSeg=zeros(HL{5},4);
            for r=1:HL{5}
                HoleSeg(Counter+r,1)=r+Counter;
                if r==HL{5}
                    HoleSeg(Counter+r,2)=HoleNd(Counter+r,1);
                    HoleSeg(Counter+r,3)=HoleNd(Counter+1,1);
                else
                    HoleSeg(Counter+r,2)=HoleNd(Counter+r,1);
                    HoleSeg(Counter+r,3)=HoleNd(Counter+r+1,1);
                end
                HoleSeg(Counter+r,4)=-1;
            end
            % Counter update
            Counter=Counter+HL{5};
            % Hole information for .poly file in Triangle
            Hole(k,1)=k;
            Hole(k,2)=Ori_cir_x;
            Hole(k,3)=Ori_cir_y;
        else
            error('Other holes than rectangle or circle is currently not available!');
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generating nodes and segments on outer rectangular boundaries
% Tuning the node density on the outer channel boundary
L=H*RatioL2H;
[N_H, N_L]=PreTri_BndAdj_rec(H,L,N_H,N_L,Fd);
% Generate the outer boundary nodes clockwise, the first node is the very first node right to the top left node. Then the top left node will be the last one 
dx=L/(N_L-1);
dy=H/(N_H-1);
N=2*N_L+2*N_H-4; % Total number of nodes on the outer boundaries.
OutNd=PreTri_Generate_OutBN(X1,X2,Y1,Y2,dx,dy,N_L,N_H,N,Fcob);
% Generate the segments on outer boundaries
OutSeg=zeros(N,4);
for r=1:N
    OutSeg(r,1)=r;
    if r==N
        OutSeg(r,2)=OutNd(r,1);
        OutSeg(r,3)=OutNd(1,1);
    else
        OutSeg(r,2)=OutNd(r,1);
        OutSeg(r,3)=OutNd(r+1,1);
    end
    OutSeg(r,4)=1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Combining the nodes and segments of inner and outer boundaries into one
% variable Nd and Seg
if Num_holes==0
    M=0;
else
    M1=length(HoleNd(:,1));
    M2=length(HoleSeg(:,1));
    if M1~=M2
        error('The number of hole segments and hole nodes does not match!');
    else
        M=M1;
    end
end
Nd=zeros(M+N,4);
Seg=zeros(M+N,4);
for r=1:M+N
    if r<=M
        Nd(r,:)=HoleNd(r,:);
        Seg(r,:)=HoleSeg(r,:);
    else
        Temp_nd=OutNd(r-M,:);
        Temp_nd(1,1)=Temp_nd(1,1)+M;
        Nd(r,:)=Temp_nd;
        Temp_seg=OutSeg(r-M,:);
        Temp_seg(1,1)=Temp_seg(1,1)+M;
        Temp_seg(1,2)=Temp_seg(1,2)+M;
        Temp_seg(1,3)=Temp_seg(1,3)+M;
        Seg(r,:)=Temp_seg;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check the combined nodes and segments
L1=length(Nd(:,1));
L2=length(Seg(:,1));
if L1~=L2
    error('The number of combined hole segments and hole nodes does not match!');
else
    L=L1;
end
%% Generate the interior nodes if interior mesh control is required.
%%%%%%%%%%%%%% Buffer layler for the inner boundaries %%%%%%%%%%%%%%%%%%%%
if Num_holes~=0
    if FBuf_inner==0
        ;
    else
        % Check whether the buffer layer will exceed the outer boundary, if
        % exceed, decrease the number of buffer layers
        if HL{1}==0
            ;
        elseif HL{1}==1 % Circular holes
            h=(dx+dy)/2; % The size in the surrounding area
            FBuf_inner=round(anylog(1+pi*sqrt(3)/HL{5},h*HL{5}/2/pi/HL{2})); % Determine the # of
            msg=['The total number of buffer layer is ', num2str(FBuf_inner), '.' ];
            disp(msg);
            % buffer layers that will smoothly merge with the outer region
            ex=1;
            while ex
                r_buff=HL{2}*(1+pi*sqrt(3)/HL{5})^FBuf_inner;
                if Ori_cir_x+r_buff>=X2 || Ori_cir_x-r_buff<=X1 || Ori_cir_y+r_buff>=Y2 || Ori_cir_y-r_buff<=Y1
                    disp('The buffer layer is beyond the outer boundaries! Decrease the number of buffer layers!');
                    if FBuf_outer==0
                        FBuf_inner=FBuf_inner-3;
                    else % The mesh also has buffer layers for the outer boundary
                        FBuf_inner=FBuf_inner-5;
                    end
                    msg=['The total number of buffer layer is ', num2str(FBuf_inner), '.' ];
                    disp(msg);
                else
                    ex=0;
                end
            end
        else
            ;
        end
        % Create buffer layers
%         if Fw==0
%             ;
%         else
%             if mod(FBuf_inner,2)~=0 % The number of buffer layer is odd
%                 FBuf_inner=FBuf_inner+1;
%                 msg=['Due to the refinement in the wake region, one more buffer layer is added, ', 'the total number of buffer layer is ', num2str(FBuf_inner), '.' ];
%                 disp(msg);
%             else
%                 ;
%             end
%         end
        for i=1:Num_holes
            if Num_holes>1
                error('Currently the buffer layers could only added to only one hole!');
            end
            
            if Num_holes==1
                HL=holes_para;
            else
                HL=holes_para{i};
            end
            
            for j=1:FBuf_inner
                if j==1
                    BASE=Nd(1:M,:); % The boundary nodes on the hole
                else
                    BASE=buff_nd_inner;
                end
                [buff_nd_inner,buff_seg_inner,S_inner]=PreTri_buff_inner(BASE,FBuf_inner,HL,Origin);
                % Merge the buffer layer with existing files
                for r=1:S_inner
                    Temp_nd=buff_nd_inner(r,:);
                    Temp_nd(1,1)=Temp_nd(1,1)+L;
                    Nd(r+L,:)=Temp_nd;
                    Temp_seg=buff_seg_inner(r,:);
                    Temp_seg(1,1)=Temp_seg(1,1)+L;
                    Temp_seg(1,2)=Temp_seg(1,2)+L;
                    Temp_seg(1,3)=Temp_seg(1,3)+L;
                    Seg(r+L,:)=Temp_seg;
                end
                %% Check the combined nodes and segments
                L1=length(Nd(:,1));
                L2=length(Seg(:,1));
                if L1~=L2
                    error('The number of combined hole segments and hole nodes does not match!');
                else
                    L=L1;
                end
            end
        end
    end
end
%%%%%%%%%%% Interior nodes in the wake region after the holes %%%%%%%%%%%%%

%%%%%%%%%%%%%% Buffer layler for the outer boundaries %%%%%%%%%%%%%%%%%%%%
if FBuf_outer==0
    ;
else
    if FBuf_inner==0
        FBuf_outer=round(min(N_H,N_L)/3);
    else % The mesh also has buffer layers fot he inner boundary
        if Num_holes==0
            ;
        elseif Num_holes==1
            if HL{1}==0 % rectangular hole
                ;
            elseif HL{1}==1 % Circular holes
                % Determine number of outer layers
                % clearence between left outer boundary and the hole
                c1=abs(X1-(Ori_cir_x-r_buff));
                % clearence between right outer boundary and the hole
                c2=abs(X2-(Ori_cir_x+r_buff));
                % clearence between bottom outer boundary and the hole
                c3=abs(Y1-(Ori_cir_y-r_buff));
                % clearence between top outer boundary and the hole
                c4=abs(Y2-(Ori_cir_y+r_buff));
                c=[c1,c2,c3,c4];
                c_min=min(c);
                if c_min<=(dx+dy)/2*1.5
                    error('The clearance between the outer boundary and the buffer layers of the inner boundary is too small! Decrease the number of buffer layers for the inner boundary!');
                end
                FBuf_outer=round(c_min/((dx+dy)/2*1.5));
                if FBuf_outer>2
                    FBuf_outer=2;
                end
            else
                ;
            end
        else
            ;
        end
    end
    % Create the buffer layers for the outer boundary
    for j=1:FBuf_outer
        if j==1
            BASE=Nd(M+1:M+N,:); % The boundary nodes on the hole
        else
            BASE=buff_nd_outer;
        end
        [buff_nd_outer,buff_seg_outer,S_outer]=PreTri_buff_outer(BASE);
        % Merge the buffer layer with existing files
        for r=1:S_outer
            Temp_nd=buff_nd_outer(r,:);
            Temp_nd(1,1)=Temp_nd(1,1)+L;
            Nd(r+L,:)=Temp_nd;
            Temp_seg=buff_seg_outer(r,:);
            Temp_seg(1,1)=Temp_seg(1,1)+L;
            Temp_seg(1,2)=Temp_seg(1,2)+L;
            Temp_seg(1,3)=Temp_seg(1,3)+L;
            Seg(r+L,:)=Temp_seg;
        end
        %% Check the combined nodes and segments
        L1=length(Nd(:,1));
        L2=length(Seg(:,1));
        if L1~=L2
            error('The number of combined hole segments and hole nodes does not match!');
        else
            L=L1;
        end
    end
end

%% Plot the combined nodes and segments
figure(1);
clf;
for r=1:L
    if r<=M
        plot(Nd(r,2),Nd(r,3),'Marker', 'o','Markersize',4, 'MarkerFaceColor', 'blue');
    elseif r<=M+N
        plot(Nd(r,2),Nd(r,3),'Marker', 'o','Markersize',4, 'MarkerFaceColor', 'red','MarkerEdgeColor', 'red');
    else
        plot(Nd(r,2),Nd(r,3),'Marker', 'o','Markersize',2, 'MarkerFaceColor', 'black','MarkerEdgeColor', 'black');
    end
    hold on;
end
axis equal tight;
msg1=['The single line for nodes is:    ', num2str(L), ' 2 ', ' 0 ', ' 1 '];
disp(msg1);
msg2=['The single line for segments is:    ', num2str(L),  ' 1 '];
disp(msg2);
msg3=['The single line for holes is:    ', num2str(Num_holes)];
disp(msg3);