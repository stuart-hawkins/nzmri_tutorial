% Compute a radiating wavefunction expansion from far field data
%
%   Note: This code requires TMATROM, which is available at http://www.romapp.org. 
%
%   pp = farfield2radiatingwavefunctionexpansion(n,kwave,[]) returns
%   observation directions (as angles) at which the far field is required
%   to compute the radiating wavefunction expansion.
%
%   a = farfield2radiatingwavefunctionexpansion(n,kwave,f) computes the
%   radiating wavefunction expansion of order n whose farfield approximates
%   f. Here f is a vector of far field values at the observation angles pp.
%
%   a/pp = farfield2radiatingwavefunctionexpansion(...,m) where m is an
%   integer uses m times as many observation angles. This can increase the
%   accuracy.
%
% Algorithm: this uses the discrete orthogonal projection given by (16) in
% Reference [1].
%
% References:
%
%   [1] Approximation of radiating waves in the near-field: error estimates
%   and application to a class of inverse problems, J. Barkhan, M. Ganesh
%   and S. C. Hawkins, J. Comput. Appl. Math. vol 401 pp 113769.
% 
% Stuart C. Hawkins - 7 May 2024

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


function varargout = farfield2radiatingwavefunctionexpansion(nmax,kwave,farfield,quadmult,docheck)

% set default for quadmult
if nargin<4 || isempty(quadmult)
    quadmult = 2;
end

%------------------------------------------
% compute quadrature points and weights
%------------------------------------------

% compute rectangle rule points
np = 2 * quadmult * (nmax+1);
pp = 2*pi*(0:np-1).'/np;
pw = 2*pi/np * ones(np,1);

%------------------------------------------
% return quadrature points if required
%------------------------------------------

if nargin<3 || isempty(farfield)
    
    % return quadrature points
    varargout{1} = pp;
    return
    
end

%------------------------------------------
% compute coefficients
%------------------------------------------

% check size of farfield matches
if size(farfield,1) ~= np || size(farfield,2) ~=1
    error('farfield must be %d x %d array',np,1);
end

% initialise coefficient array
cof = zeros(2*nmax+1,1);

% compute the coefficients
for j=-nmax:nmax
    
    cof(j+nmax+1) = sum(pw.*exp(-1i*j*pp).*farfield);
    
end

% apply the scaling
s = 0.25*(1+1i)*sqrt(kwave/pi)*1i.^abs(-nmax:nmax).';

cof = s.*cof;

% compute radiating wavefunction expansion
varargout{1} = radiatingwavefunctionexpansion(nmax,0,kwave,cof);

%------------------------------------------
% check the error in the coefficients
%------------------------------------------

if nargin>4

    % compute the far field
    check = varargout{1}.evaluateFarField(exp(1i*pp));

    % compute and return the error
    varargout{2} = sqrt(sum(2*pi*abs(check-farfield).^2)/length(pp));
    
end
