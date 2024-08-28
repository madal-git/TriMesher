function [n,star]=PostTri_findstar(nd,ele)
% function [n,star]=find_star(node,ele) find the triangles that form a star
% structure for the given node
% ele---The elelment file from Triangle
% nd---The node order number of given node
% n---The total number of triangles in each star structure
% star---The row vector containing the order number of all triangles in
% each star structure

n=0;

M=length(ele(:,1));
for r=1:M
    E=ele(r,:);
    if E(1,2)==nd || E(1,3)==nd || E(1,4)==nd
        n=n+1;
        star(1,n)=r;
    end
end
% check whether the found duplicate triangles in the star structure
for t=1:n
    exclude=star(1,t);
    star_compare=zeros(1,n-1);
    a=0;
    for s=1:n
        if star(1,s)~=exclude
            a=a+1;
            star_compare(1,a)=star(1,s);
        end
    end
    if a~=n-1
        error('The star structure for comparison is not correct!');
    end
    for s=1:a
        if exclude==star_compare(1,s)
            error('Found duplicate triangles in the star structure');
        end
    end
end