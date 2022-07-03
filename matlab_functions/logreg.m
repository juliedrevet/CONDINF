function [b] = logreg(fun,varargin)
%  LOGREG  Perform binary logistic regression
%  based on Valentin Wyart's function - 09/2015
%  
%  Usage: [b] = LOGREG(fun,x,y,w,[x],[y],[w])
%
%  where fun   - psychometric curve equation
%        x     - design matrix (input)
%        y     - binary responses (output)
%        w     - weight attributed to x,y relation
%
%  (x,y,w) works as a triplet, function can take several triplets as argument
%
%  The function returns maximum-likelihood parameter estimates b obtained via
%  gradient descent using fmincon.


% check input arguments
if (nargin < 4) || (mod(numel(varargin),3)) % 
    error('Missing input arguments!');
end

sizes = zeros(1,numel(varargin)/3);
for k = 1:(numel(varargin)/3)
    X{k} = varargin{(k-1)*3+1};
    Y{k} = varargin{(k-1)*3+2};
    W{k} = varargin{(k-1)*3+3};  
    
    sizes(k) = length(X{k});
    
    if any(size(X{k}) ~= size(Y{k}))
        error('Mismatching input/output arrays!');
    end
    if ~all(size(W{k}) == 1)
        if any(size(W{k}) ~= size(Y{k}))
            error('Mismatching input/output arrays!');
        end
    end
end

x = [];
y = [];
w = [];

for k = 1:(numel(varargin)/3)
    x = [x; reshape(X{k},sizes(k),1)];
    y = [y; reshape(Y{k},sizes(k),1)];
    if all(size(W{k}) == 1)
        w = [w; W{k}*ones(sizes(k),1)];
    else
        w = [w; reshape(W{k},sizes(k),1)];
    end
end




% ensure binary output array
y = logical(y);
% y = reshape(y,numel(y),1);
% x = reshape(x,numel(x),1);

% get maximum-likelihood parameter estimates
n = size(x,2);

% fun = @(p,x) normcdf(p(:,1)'+p(:,2)'.*x)*(1-p(:,3)') + p(:,3)'/2;
fun = @(pp,xx) pp(:,3)/2+normcdf(pp(:,1)+pp(:,2)*xx)*(1-pp(:,3));
ub = repmat([0 1000 1],n,1);
lb = repmat([-100 0 0],n,1);
X0 = [-10 20 0.5];


Eq = @(pp)-get_llh(pp);

b = fmincon(Eq,repmat(X0,n,1),[],[],[],[],lb,ub,[],...
 optimset('Display','notify','FunValCheck','on','Algorithm','interior-point','TolX',1e-20,'MaxFunEvals',1e6)); % if commented, warning!

    function [llh] = get_llh(p)
        % get logistic regression model log-likelihood

        plh = fun(p,x);
        plh(y==0) = 1-plh(y==0);
        llh = sum(w.*log(max(plh,eps)));
    end

end