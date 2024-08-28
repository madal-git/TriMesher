function [X,Y,N]=node_combine(X_B,Y_B,N_B,X_I_B,Y_I_B,N_I_B,X_I_AB,Y_I_AB,N_I_AB)
% Put the coordinates of all nodes into X and Y.
% X_B,Y_B and N_B are defined in the function B_node.
% X_I_B,Y_I_B and N_I_B are defined in the function I_nodeB.
% X_I_AB,Y_I_AB and N_I_AB are defined in the function I_nodeAB.
% X and Y are passed to main function.
% N is the total number of nodes.

N=N_B+N_I_B+N_I_AB;
X=zeros(N,1);
Y=zeros(N,1);

for i=1:N;
    if i<=N_I_B;
        X(i)=X_I_B(i); % Matlab naturely proceed in column in 2D matrix
        Y(i)=Y_I_B(i); % Matlab naturely proceed in column in 2D matrix
    elseif i<=N_I_B+N_I_AB;
        X(i)=X_I_AB(i-N_I_B); % Matlab naturely proceed in column in 2D matrix
        Y(i)=Y_I_AB(i-N_I_B); % Matlab naturely proceed in column in 2D matrix
    else
        X(i)=X_B(i-N_I_B-N_I_AB);
        Y(i)=Y_B(i-N_I_B-N_I_AB);
    end;
end;