function f=mynormpdf(z,mu,sig)
% f=mynormpdf(z,mu,sig)
% INPUT: z = values to estimate PDF at
%          mu = mean value
%          sig = standard deviation
% OUTPUT: f = probability density function for a normal/Gaussian distribution

A=1/(sig*sqrt(2*pi)); % constant to make the PDF integrate to 1 over +/−inf
B=(z - mu).^2; % squared deviation from mean to give nomerator in exponent
C=2*sig.^2; % 2 times the standard deviation squared to give denominator in exponent
f=A*exp(-B./C); % normal PDF
