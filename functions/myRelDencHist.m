function [N,binx] = myRelDencHist(D,nbins)
% [N,binx] = myRelDencHist(D,nbins)
% INPUT: 
%       D = Dataset
%       nbins = number of bins
% OUTPUT:
%       N = normalized histogram bin counts
%       binx = bin definitions for the relative dencity histogram

[N,binx] = hist(D,nbins); %making a histogram of D with 'nbins' number of bins
dx = binx(2)-binx(1); %calculating the width of bins
Area = sum(N*dx); %calculating the area under the historgram
N = N/Area; %normalizing the number of counts in each bins
