function [X_I_B,Y_I_B,N_I_B]=I_nodeB(N_H,N_L,X1,Y1,h)
% function [X_I_B,Y_I_B,N_I_B]=I_nodeB(N_H,N_L,X1,Y1,h)
% generates the interior nodes that are aligned with boundary nodes.
% N_H is defined in the main function.
% N_L is defined in the function B_node.
% X1 and Y1 are defined in function ini_tri.
% h is defined in function B_node.
% X_I_B is the x coordinates of interior nodes that are aligned with
% boundary nodes.
% Y_I_B is the y coordinates of interior nodes that are aligned with
% boundary nodes.
% N_I_B is the number of interior nodes that are aligned with
% boundary nodes.

N_I_B=0;
for i=1:N_L-2;
    for j=1:N_H-2;
        X_I_B(i,j)=X1+i*h;
        Y_I_B(i,j)=Y1+j*h;
        N_I_B=N_I_B+1;
%         Message=['Node ' , num2str(N_I_B) , ' is created in interior region.'];
%         disp(Message);
    end;
end;