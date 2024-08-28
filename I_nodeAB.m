function [X_I_AB,Y_I_AB,N_I_AB]=I_nodeAB(N_H,N_L,X1,Y1,h,N_I_B)
% function [X_I_AB,Y_I_AB,N_I_AB]=I_nodeAB(N_H,N_L,X1,Y1,h)
% generates the interior nodes that are not aligned with boundary nodes.
% N_H is defined in the main function.
% N_L is defined in the function B_node.
% X1 and Y1 are defined in function ini_tri.
% h is defined in function B_node.
% N_I_B is defined in the function I_nodeB.
% X_I_AB is the x coordinates of interior nodes that are not aligned with
% boundary nodes.
% Y_I_AB is the y coordinates of interior nodes that are not aligned with
% boundary nodes.
% N_I_AB is the number of interior nodes that are not aligned with
% boundary nodes.

N_I_AB=0;
for i=1:N_L-1;
    for j=1:N_H-1;
        X_I_AB(i,j)=X1+h/2+(i-1)*h;
        Y_I_AB(i,j)=Y1+h/2+(j-1)*h;
        N_I_AB=N_I_AB+1;
%         Message=['Node ' , num2str(N_I_B+N_I_AB) , ' is created in interior region'];
%         disp(Message);
    end;
end;