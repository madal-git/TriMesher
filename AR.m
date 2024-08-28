function ratio=AR(P1,P2,P3)
%P1, P2 and P3 are the coordinates of the triangle vertices
%Each Input should be placed in a columm vector
%Output ratio will be single numeric value of aspect ratio for triangle

L1 = dis(P1,P2); %Distance from point one to point two
L2 = dis(P1,P3); %Distance from point one to point three
L3 = dis(P2,P3); %Distance from point two to point three

lengthmax = max(max(L1,L2),L3);
lengthmin = min(min(L1,L2),L3);

ratio = lengthmax/lengthmin;


