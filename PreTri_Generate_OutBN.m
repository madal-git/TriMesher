function OutNd = PreTri_Generate_OutBN(X1,X2,Y1,Y2,dx,dy,N_L,N_H,N,Fcob)

% function OutNd = PreTri_Generate_OutBn(X1,X2,Y1,Y2,dx,dy,N_L,N_H,N,Fcob)
% generates boundary nodes on the outer boundaries of the computational
% domain
% Fcob is the flag for the shape of outee boundaries. 0---Rectanglar; non-zero integers---other customized shape

OutNd=zeros(N,4);

if Fcob==0 % Rectangular
    for r=1:N
        OutNd(r,1)=r;
        if r<=N_L-2
            OutNd(r,2)=X1+r*dx;
            OutNd(r,3)=Y2;
        elseif r==N_L-1
            OutNd(r,2)=X2;
            OutNd(r,3)=Y2;
        elseif r<=N_L-1+N_H-2
            OutNd(r,2)=X2;
            OutNd(r,3)=Y2-(r-N_L+1)*dy;
        elseif r==N_L-1+N_H-1
            OutNd(r,2)=X2;
            OutNd(r,3)=Y1;
        elseif r<=N_L-1+N_H-1+N_L-2
            OutNd(r,2)=X2-(r-N_L+1-N_H+1)*dx;
            OutNd(r,3)=Y1;
        elseif r==N_L-1+N_H-1+N_L-1
            OutNd(r,2)=X1;
            OutNd(r,3)=Y1;
        elseif r<N
            OutNd(r,2)=X1;
            OutNd(r,3)=Y1+(r-N_L+1-N_H+1-N_L+1)*dy;
        else
            OutNd(r,2)=X1;
            OutNd(r,3)=Y2;
        end
        OutNd(r,4)=1;
        %     figure(1);
        %     plot(OutNd(r,2),OutNd(r,3),'Marker', 'o','Markersize',4, 'MarkerFaceColor', 'red','MarkerEdgeColor', 'red');
        %     axis equal tight;
        %     hold on;
    end
elseif Fcob==1 % A half circular hump in the middle of the bottom surface
    n=3;
    humP_radius=n*dx;
    %% Initial rectangular shape
    for r=1:N
        OutNd(r,1)=r;
        if r<=N_L-2 % Top surface
            OutNd(r,2)=X1+r*dx;
            OutNd(r,3)=Y2;
        elseif r==N_L-1
            OutNd(r,2)=X2;
            OutNd(r,3)=Y2;
        elseif r<=N_L-1+N_H-2 % Right surface
            OutNd(r,2)=X2;
            OutNd(r,3)=Y2-(r-N_L+1)*dy;
        elseif r==N_L-1+N_H-1
            OutNd(r,2)=X2;
            OutNd(r,3)=Y1;
        elseif r<=N_L-1+N_H-1+N_L-2
            OutNd(r,2)=X2-(r-N_L+1-N_H+1)*dx; % Bottom surface
            OutNd(r,3)=Y1;
        elseif r==N_L-1+N_H-1+N_L-1
            OutNd(r,2)=X1;
            OutNd(r,3)=Y1;
        elseif r<N % Left surface
            OutNd(r,2)=X1;
            OutNd(r,3)=Y1+(r-N_L+1-N_H+1-N_L+1)*dy;
        else
            OutNd(r,2)=X1;
            OutNd(r,3)=Y2;
        end
        OutNd(r,4)=1;
    end

    %% Find the node that is closet to the middle point
    Dis_2_mid=zeros(N,1);
    mid_coord=[(X2-X1)/2;Y1];
    for r=1:N
        if r<=N_L-1+N_H-1+N_L-2 && r>N_L-1+N_H-1
            Dis_2_mid(r)=dis(mid_coord,[OutNd(r,2);OutNd(r,3)]);
        end
    end
    Dis_2_mid_min=min(Dis_2_mid(N_L-1+N_H:N_L-1+N_H-1+N_L-2));
    for r=1:N
        if r<=N_L-1+N_H-1+N_L-2 && r>N_L-1+N_H-1
            if single(dis(mid_coord,[OutNd(r,2);OutNd(r,3)]))==single(Dis_2_mid_min)
                break;
            end
        end
    end
    Circle_center_coord=[OutNd(r,2);OutNd(r,3)];
    Influence_range=r-(n-1):1:r+(n-1);
    Num_of_pie_slices=2*n;
    Angle_each_pie_slice=pi/Num_of_pie_slices;
    %% Create new coordinates for the nodes on the half-cricle hump
    counter=0;
    for r=1:N
        if r<=Influence_range(end) && r>=Influence_range(1)
            counter=counter+1;
            OutNd(r,2)=Circle_center_coord(1,1)+humP_radius*cos(counter*Angle_each_pie_slice);
            OutNd(r,3)=Circle_center_coord(2,1)+humP_radius*sin(counter*Angle_each_pie_slice);
        end
    end
    
else
    error('Wrong flag for the outer boundary shape!');
end