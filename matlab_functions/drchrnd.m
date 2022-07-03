function [r] = drchrnd(a,n)
%  DRCHRND  Random arrays from the Dirichlet distribution
a = reshape(a,1,[]);
p = length(a);
r = gamrnd(repmat(a,n,1),1,n,p);
r = r./repmat(sum(r,2),1,p);
end