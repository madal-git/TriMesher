function skewness = SK(P1,P2,P3)
%P1, P2 and P3 are the coordinates of the triangle vertices
%Each Input should be placed in a columm vector
%Output ratio will be single numeric value of aspect ratio for triangle

u1 = P1-P2;
u2 = P1-P3;
u3 = P2-P3;

theta1 = angle(u1,u2);
theta2 = angle(u1,u3);
theta3 = angle(u2,u3);

thetamax = max(max(theta1,theta2),theta3);
thetamin = min(min(theta1,theta2),theta3);
skewness = max(((thetamax-60)/120),((60-thetamin)/60));