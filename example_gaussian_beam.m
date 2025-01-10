% Example script for fast acoustic sound-soft multiple-scattering simulation
%
%   example_gaussian_beam computes the scattered field for the point-source
%   problem (39) in [1], and scattering by a plane wave or Gaussian beam.
%   The exact solution is known for the point-source problem, and so the
%   error is computed and printed in the command window. The scattered
%   field is computed using fast_multiple_scattering_soft.m, which
%   implements the algorithm in [1].
%
% The scatterer positions, wavenumber etc are set in
% fast_multiple_scattering_hanmer_springs.m
%
% Algorithm: this code uses the fast Stage 3 algorithm in [1] and the
% TMATROM package from http://www.romapp.org.
%
% See also: example_sound_soft, example_transmission, plane_wave.
%
% References:
%
% [1] A fast algorithm for the two-dimensional Helmholtz transmission
% problem with large multiple scattering configurations, S. C. Hawkins and
% M. Ganesh, J. Acoust. Soc. Am 2024.
%
% Stuart C. Hawkins - 10 Jan 2025

% Copyright 2025 Stuart C. Hawkins
% 	
% This file is part of nzmri_tutorial.
% 
% nzmri_tutorial is free software: you can redistribute it
% and/or modify	it under the terms of the GNU General Public License as
% published by the Free Software Foundation, either version 3 of the License,
% or (at your option) any later version.
% 
% nzmri_tutorial is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with nzmri_tutorial. If not, see <http://www.gnu.org/licenses/>.


clear all
close all

% set type of incident wave... options are POINTSOURCE, PLANEWAVE or
% GAUSSIANBEAM
type = 'GAUSSIANBEAM';

% create a grid of points at which to compute the total field
x = linspace(-5,15,100);
y = linspace(-20,20,200);
[xx,yy] = meshgrid(x,y);
zz = xx + 1i*yy;

% call the function to compute the far field
[u,I,target] = fast_multiple_scattering_hanmer_springs(zz,type);

%%

% plot the exterior field
figure
surf(xx,yy,zeros(size(u)),real(u))
shading interp
axis equal
view([0 90])
hold on
plot(real(target),imag(target),'kx','markersize',10,'linewidth',4)
hold off

% print out the intensity at the target
fprintf('Intensity at the target is |u|^2 = %0.3e\n',I);