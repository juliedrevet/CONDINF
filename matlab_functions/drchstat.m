function [m,v] = drchstat(a)
%  DRCHSTAT  Means and variances for the Dirichlet distribution
A = sum(a);
m = a./A;
v = a.*(A-a)./((A+1)*A^2);
end