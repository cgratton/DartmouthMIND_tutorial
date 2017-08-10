function [matrix r kden] = threshold_the_matrix(matrix,threshold)
% based on a script originally from JDP
% density threshold a matrix
% CG edits

d=size(matrix);
numpossibleedges=d(1)*(d(1)-1)/2;

edgesleft=ceil(threshold*numpossibleedges);
matrix=triu(matrix,1);
[v i]=sort(matrix(:),'descend');
keptedges=zeros(numpossibleedges,1);
keptedges(1:edgesleft)=1;
keptedges=logical(keptedges);
v(~keptedges)=0;

matrix(i)=v;
matrix=reshape(matrix,d);
matrix=max(matrix,matrix');
r=v(edgesleft);
kden=edgesleft/numpossibleedges;

end