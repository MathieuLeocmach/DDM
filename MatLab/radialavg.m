function [Zr, R] = radialavg(z,m)
% RADIALAVG	 Radially averaqe 2D square matrix Z into M bins
%
% [Zr, R] = radialavg(z,m)
%
% [Zr, R] = RADIALAVG(z,m) computes the average along the radius of a unit
% circle inscribed in the square matrix z. The average is computed in M
% bins. The radial average is not computed beyond the unit circle, in the
% corners of the matrix z. The radial average is returned in Zr and the
% mid-points of the M bins are returned in vector R.
%
% Example
%	N=101;
%	[X,Y] = meshgrid(-1:2/(N-1):1);
%	z = 1-sqrt(X.^2 + Y.^2);
%	m=(N-1)/2+1;
%	[Zr,R] = radialavg(z,m);
%	figure;plot(R,Zr,'.-');
%
% INPUT
% z = square input matrix to be radially averaged
% m = number of bins in which to compute radial average
%
% OUTPUT
% Zr = radial average of length m
% R  = m locations of Zr (i.e. midpoints of the m bins)
% 
% See also linspace, meshgrid

% (c) 2014 David J. Fischer | fischer@shoutingman.com
% 4/4/14 DJF first working version
% 5/2/14 DJF documentation & radialavg_tester.m to simply demonstrate use

% radial distances r over grid of z
N = size(z,1);
[X,Y] = meshgrid(-1:2/(N-1):1);
r = sqrt(X.^2+Y.^2);

% equi-spaced points along radius which bound the bins to averaging radial values
% bins are set so 0 (zero) is the midpoint of the first bin and 1 is the last bin
dr = 1/(m-1);
rbins = linspace(-dr/2,1+dr/2,m+1);

% radial positions are midpoints of the bins
R =(rbins(1:end-1)+rbins(2:end))/2;

Zr = zeros(1,m); % vector for radial average

% loop over the bins, except the final (r=1) position
for j=1:m-1
	% find all matrix locations whose radial distance is in the jth bin
	bins = r>=rbins(j) & r<rbins(j+1);
	
	% count the number of those locations
	n = sum(sum(bins));
	if n~=0
		% average the values at those binned locations
		Zr(j) = sum(z(bins))/n;
	else
		% special case for no bins (divide-by-zero)
		Zr(j) = NaN;
	end
end

% special case the last bin location to not average Z values for
% radial distances in the corners, beyond R=1
bins = r>=rbins(m) & r<=1;
n = sum(sum(bins));
if n~=0
	Zr(m) = sum(z(bins))/n;
else
	Zr(m) = NaN;
end
