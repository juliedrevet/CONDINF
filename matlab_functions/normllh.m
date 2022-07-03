function y = normllh(x,mu,sigma)
%  NORMLLH  Normal log-likelihood
%
%  This is a copy of NORMPDF returning the logarithm of the probability density
%  function to avoid numerical precision issues (nulling low probabilities).
%
%  Valentin Wyart <valentin.wyart@ens.fr>

if nargin<1
    error(message('stats:normpdf:TooFewInputs'));
end
if nargin < 2
    mu = 0;
end
if nargin < 3
    sigma = 1;
end

% Return NaN for out of range parameters.
sigma(sigma <= 0) = NaN;

try
    y = -0.5 * ((x - mu)./sigma).^2 - log(2*pi)/2 - log(sigma);
catch
    error(message('stats:normpdf:InputSizeMismatch'));
end
