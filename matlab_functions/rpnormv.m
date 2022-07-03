function [x,c] = rpnormv(m,s)
%  RPNORMV  Sample random numbers from positive normal distributions
%
%  Usage: [x,c] = RPNORMV(m,s)
%
%  where m is an array of means for the normal distributions and s an array of
%  standard deviations for the normal distributions. The function has been
%  vectorized for maximal efficiency. The output sample array x will match in
%  size the larger of the two input arrays. The second output array c is the
%  number of rejected propositions for each sample.
%
%  This function is a vectorized version of the rpnorm function which draws
%  multiple samples from a single positive normal distribution. RPNORMV allows
%  to draw single samples from multiple positive normal distributions with no
%  loss in speed.
%
%  Reference:
%      V. Mazet, D. Brie, J. Idier, "Simulation of Positive Normal Variables
%      using several Proposal Distributions", IEEE Workshop on Statistical
%      Signal Processing 2005, July 17-20 2005, Bordeaux, France
%
%  Valentin Wyart <valentin.wyart@ens.fr>, 10/2015

% check input arguments
if nargin < 2
    s = 1;
end
if nargin < 1
    error('Missing input arguments!');
end
if ~all(size(m) == size(s)) && ~(isscalar(m) || isscalar(s))
    error('Mismatching input arguments!');
end
if any(isinf(m) | isnan(m)) || any(isinf(s) | isnan(s) | s <= 0)
    error('Invalid input arguments!');
end

% store input size
if ~isscalar(m)
    siz = size(m);
else
    siz = size(s);
end

% reshape input arrays
if isscalar(m)
    m = repmat(m,siz);
elseif isscalar(s)
    s = repmat(s,siz);
end

% columnize input arrays
m = m(:);
s = s(:);

% initialize output array
N = prod(siz);
x = nan(N,1);
c = zeros(N,1);

% intersections
A  = 1.136717791056118;
mA = (1-A^2)/A.*s;
mC = s.*sqrt(pi/2);

while any(isnan(x))
    
    z = nan(N,1);
    rho = nan(N,1);
    
    % ignore already accepted samples
    rho(~isnan(x)) = -1;
    
    % 4/ exponential distribution
    i = isnan(rho) & m < mA; n = nnz(i);
    if n > 0
        a = (-m(i)+sqrt(m(i).^2+4*s(i).^2))/2./s(i).^2;
        z(i) = -log(1-rand(n,1))./a;
        rho(i) = exp(-(z(i)-m(i)).^2/2./s(i).^2-a.*(m(i)-z(i)+a.*s(i).^2/2));
    end
    
    % 3/ normal distribution truncated at the mean
    i = isnan(rho) & m <= 0; n = nnz(i);
    if n > 0
        z(i) = abs(randn(n,1)).*s(i)+m(i);
        rho(i) = (z(i) >= 0);
    end
    
    % 2/ normal distribution coupled with uniform
    i = isnan(rho) & m < mC; n = nnz(i);
    if n > 0
        r = (rand(n,1) < m(i)./(m(i)+sqrt(pi/2)*s(i)));
        u = rand(n,1).*m(i);
        g = abs(randn(n,1).*s(i))+m(i);
        z(i) = r.*u+(1-r).*g;
        rho(i) = r.*exp(-(z(i)-m(i)).^2/2./s(i).^2)+(1-r);
    end
    
    % 1/ normal distribution
    i = isnan(rho); n = nnz(i);
    if n > 0
        z(i) = randn(n,1).*s(i)+m(i);
        rho(i) = (z(i) >= 0);
    end
    
    % accept valid propositions
    i = (rand(N,1) <= rho); % accepted propositions
    j = ~i & isnan(x); % rejected propositions
    x(i) = z(i);
    c(j) = c(j)+1;
    
end

% reshape output array
x = reshape(x,siz);
c = reshape(c,siz);

end