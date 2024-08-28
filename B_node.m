function [X_B,Y_B,N_L,N_B,h]=B_node(N_H,L,RatioLH,X1,X2,Y1,Y2)
% [X_B,Y_B,N_B]=B_node(N_H,L,RatioLH) generates the nodes
% on boundaries of the rectangular.
% N_H, L and RatioLH are defined in the main function.
% X1 is defined in the function ini_tri.
% X2 is defined in the function ini_tri.
% Y1 is defined in the function ini_tri.
% Y2 is defined in the function ini_tri.
% X_B is the x coordinates of the nodes on the boundaries.
% Y_B is the y coordinates of the nodes on the boundaries.
% N_L is the number of nodes on the length edge.
% N_B is the number of nodes on the boundaries.
% h is the spacing between two nodes.

N_L=(N_H-1)*RatioLH+1; % The number of nodes on the length edge 
h=L/(N_L-1);  % The spacing in x and y is the same.
N_B=2*N_H+2*(N_L-2); % Total number of nodes on the edges of rectangle

X_B=zeros(N_B,1);
Y_B=zeros(N_B,1);

for i=1:N_B;
    if i<=N_L;
        X_B(i)=X1+(i-1)*h;
        Y_B(i)=Y2;
    elseif i<=N_L+N_H-2;
        X_B(i)=X2;
        Y_B(i)=Y2-(i-N_L)*h;
    elseif i<=N_L+N_H+(N_L-2); 
        X_B(i)=X2-(i-(N_L+N_H-2)-1)*h;
        Y_B(i)=Y1;
    else
        X_B(i)=X1;
        Y_B(i)=Y1+(i-(2*N_L+N_H-2))*h;
    end;
end;

for i=1:N_B;
    if i<N_B;
        XX(i)=X_B(i+1);
        YY(i)=Y_B(i+1);
    else
        XX(i)=X_B(1);
        YY(i)=Y_B(1);
    end;
end;

X_B=XX;
Y_B=YY;