% Fast sound-soft multiple-scattering simulation
%
%   u = fast_multiple_scattering_hanmer_springs(z,'PLANEWAVE') computes the
%   total field u at the points z for multiple scattering of the
%   plane wave plane_wave(0,kwave) where kwave is the wavenumber (set
%   below). The preconditioned GMRES convergence history is plotted.
%
%   u = fast_multiple_scattering_hanmer_springs(z,'GAUSSBEAM') computes the 
%   total field u at the points z for multiple scattering of a
%   Gaussian beam gaussian_beam(pi/4,0,3,kwave) where kwave is the 
%   wavenumber (set below), pi/4 is the direction, 0 the focus point and 3
%   the waist radius. The preconditioned GMRES convergence history is 
%   plotted.
%
%   u = fast_multiple_scattering_hanmer_springs(z,'POINTSOURCE') computes
%   total field u for the point-source problem (39) in [1] (ignore the 
%   interior field). The exact far field is known for this problem, and so 
%   the error in the far field computed and printed in the command window.
%
%   [u,I0,target] = fast_multiple_scattering_hanmer_springs(...) computes
%   also the intensity |u|^2 at the point target. The value of target is
%   set in the code below.
%
% Note: position vectors are stored in complex format ie x+1i*y represents
% the point (x,y).
%
% Algorithm: this code uses the fast Stage 3 algorithm in [1] and the
% TMATROM package from http://www.romapp.org.
%
% See also: fast_multiple_scattering_transmission, plane_wave.
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


function [utot,I0,target] = fast_multiple_scattering_hanmer_springs(points,type)

%==================================================================
% SETUP
%==================================================================

%-------------------------------------------- 
% process parameters
%-------------------------------------------- 

% set the default for type
if nargin<2
    type = 'POINTSOURCE';
end

%-------------------------------------------- 
% set main parameters
%-------------------------------------------- 

% set the target point
target = 10i;

% shape of the scatterer - circle with radius 1
geom = obstacleCircle(1);

% wavenumber
kwave = 2;

% Nystrom discretisation parameter
n = 20;

%-------------------------------------------- 
% set the geometry of the scatterers
%-------------------------------------------- 

% number of scatterers
N = 20;

% we will put the scatterers on a straight line between x0 and x1
x0 = 10-20i;
x1 = 10+20i;

% compute the coordinates of the scatterers
pos = x0+(x1-x0)*(0:N-1)/(N-1);
       
%-------------------------------------------- 
% set dependent parameters
%-------------------------------------------- 

% incident field
switch type
    
    case 'POINTSOURCE'
        
        % point source located inside the first scatterer
        inc = point_source(pos(1),kwave);

    case 'PLANEWAVE'
        
        % plane wave with direction (1,0)
        inc = plane_wave(pi/4,kwave);

    case 'GAUSSIANBEAM'

        % Gaussian beam with direction (1,1)/sqrt(2), waist radius 6, and 
        % focus at (0,0)
        inc = gaussian_beam(pi/4,0+0i,3,kwave);
        
    otherwise
        
        error('type must be PLANEWAVE or GAUSSIANBEAM or POINTSOURCE (default).')
        
end

% setup solver to solve single-scattering problems
slvr = solverNystrom(kwave,[],geom);

% setup the solver... this creates discretisation matrix etc
slvr.setup(n);

% get the suggested order for the wavefunction expansions, based on the
% scatterer radius and wavenumber
nmax = suggestedorder(kwave,slvr.getRadius());

% get quadrature points for the far field projection (25) in [1]
ppfar = farfield2radiatingwavefunctionexpansion(nmax,kwave,[],1);

%-------------------------------------------- 
% visualise the configuration
%-------------------------------------------- 

% plot the scattering configuration
figure
plot(real(pos),imag(pos),'kx')
hold on
t = linspace(0,2*pi,50);
for j=1:N
    plot(real(pos(j))+cos(t),imag(pos(j))+sin(t),'k--')
end
hold off
axis equal

%-------------------------------------------- 
% quick check that the scatterers don't intersect
%-------------------------------------------- 

% Note: there are O(log(N)) ways to do this but this will do for small N
for j=1:N    

    % find the distance of the closest centre to the centre of scatterer j
    d = min(abs(pos(j)-pos([1:j-1,j+1:N])));

    % if the distance is less than double the radius then the scatterers
    % could intersect and we need to raise an error
    if d<2*slvr.getRadius()
        error('Two scatterers may intersect')
    end

end

%-------------------------------------------- 
% factorise the discretisation matrices... the
% factors are used for the block-Jacobi preconditioner
%-------------------------------------------- 

% get LU-factorisation of the discretisation matrix... we use this in
% the preconditioner
[L,U,P] = lu(slvr.matrix);

%==================================================================
% SOLVE THE LINEAR SYSTEM
%==================================================================

%-------------------------------------------- 
% compute the RHS
%-------------------------------------------- 

% get Nystrom points on reference domain [0,2*pi]
pp=pi*(0:2*n-1)/n;
pp=pp(:);

% apply mapping to quadrature points to get points on the scatterer... we
% use these later
% Note: these assume the scatterer is located at (0,0)
[x,y,qx,qy] = geom.geom(pp);

% initialise array to hold the right hand side
b = zeros(2*n,N);

for j=1:N

    % get points on scatterer j in in complex format for the incident field
    % evaluation
    qz = pos(j) + qx + 1i*qy;

    % evaluate the incident field and store the values
    if strcmp(type,'POINTSOURCE')
        % for point source we swap the sign because we are interested in
        % the radiation of a field induced by an interior point source...
        % see (40) in [1]
        b(:,j) = 2*inc.evaluate(qz);
    else
        b(:,j) = -2*inc.evaluate(qz);
    end
    
end

% reshape the right hand side into a vector
b = b(:);

%-------------------------------------------- 
% solve the system
%-------------------------------------------- 

% work out how many iterations to compute
itns = min(100,ceil(length(b)/2));

% solve the system using GMRES... indented functions below implement the
% matrix-vector products with the matrix and the preconditioner
[x,flag,relres,iter,resvec] = gmres(@matrix_product,b,itns,1e-8,1,...
    @preconditioner);

% plot the residual history in a figure
figure
semilogy(0:length(resvec)-1,resvec)
xlabel('iteration')
ylabel('residual norm')
title('GMRES residual history')

%==================================================================
% postprocessing
%==================================================================

%-------------------------------------------- 
% compute the total field
%-------------------------------------------- 

% reshape x so that x(:,j) are the coefficients on scatterer j
x = reshape(x,[],N);

% compute the incident field
utot = inc.evaluate(points);

% compute the far field from each scatterer in turn
for j=1:N
    
    % we use the slvr to compute the far field... we manually set the
    % coefficients inside slvr
    slvr.cof = x(:,j);

    % compute the exterior field... slvr thinks the scatterer is at the
    % origin so we adjust the points
    utot = utot + slvr.getField(points-pos(j));
    
end

% mask out the scatterers
for j=1:N
    ii = abs(points-pos(j))<slvr.getRadius();
    utot(ii) = NaN;
end

%-------------------------------------------- 
% compute the intensity at the target
%-------------------------------------------- 

% reshape x so that x(:,j) are the coefficients on scatterer j
x = reshape(x,[],N);

% compute the incident field
u0 = inc.evaluate(target);

% compute the far field from each scatterer in turn
for j=1:N
    
    % we use the slvr to compute the far field... we manually set the
    % coefficients inside slvr
    slvr.cof = x(:,j);

    % compute the exterior field... slvr thinks the scatterer is at the
    % origin so we adjust the points
    u0 = u0 + slvr.getField(target-pos(j));
    
end

% compute the intensity
I0 = abs(u0).^2;

%-------------------------------------------- 
% check the error
%-------------------------------------------- 

% for the POINTSOURCE problem the true far field is the far field from the
% point source itself

if strcmp(type,'POINTSOURCE')    

    % setup some angles at which to compute the far field
    tp = linspace(0,2*pi,1000);
    tp = tp(:);

    % - - - - - - - - - - - - - - - - - - - - 
    % compute the far field
    % - - - - - - - - - - - - - - - - - - - - 

    % reshape x so that x(:,j) are the coefficients on scatterer j
    x = reshape(x,[],N);

    % initialise array to hold the far field
    ff = zeros(size(tp,1),1);

    % compute the far field from each scatterer in turn
    for j=1:N

        % we use the slvr to compute the far field... we manually set the
        % coefficients inside slvr
        slvr.cof = x(:,j);
        tmp = slvr.getFarField(tp);

        % slvr thinks the scatterer is at (0,0) so we manually adjust the phase
        % of the far field to get the far field for a scatterer at pos(j)
        phase = exp(-1i*kwave*( cos(tp)*real(pos(j)) + sin(tp)*imag(pos(j)) ));

        % add on the far field
        ff = ff + phase.*tmp;

    end

    % - - - - - - - - - - - - - - - - - - - - 
    % calculate the reference solution... which is the field from the point
    % source
    % - - - - - - - - - - - - - - - - - - - -
    ref = inc.evaluateFarField(exp(1i*tp));

    % compute the error
    fprintf('Error for point source is %0.2e\n',max(abs(ref-ff)));
    
end

%==================================================================
% indented functions used within GMRES
%==================================================================

    %-------------------------------------------- 
    % matrix vector product using the Stage 3
    % algorithm from [1]
    %-------------------------------------------- 

    function y = matrix_product(x)

        % reshape x so that x(:,j) are the coefficients from scatterer j
        x = reshape(x,[],N);

        % initialise y, the product
        y = zeros(size(x));

        % - - - - - - - - - - - - - - - - - - - - - - - - - -
        % multiply by the block-diagonal part in (33) in [1]
        % - - - - - - - - - - - - - - - - - - - - - - - - - -
    
        % loop through all the scatterers
        for j=1:N
            % multiply by the discretisation of the block-diagonal operator
            y(:,j) = slvr.matrix * x(:,j);
        end

        % - - - - - - - - - - - - - - - - - - - - - - - - - -
        % precompute the expansion coefficients in (34) of [1]
        % - - - - - - - - - - - - - - - - - - - - - - - - - -
    
        % loop through all the scatterers
        for j=1:N
        
            % we use the slvr to compute the far field... we manually set the
            % coefficients inside slvr
            slvr.cof = x(:,j);
            ff = slvr.getFarField(ppfar);
            
            % now we get the coefficients from the far field
            rad{j} = farfield2radiatingwavefunctionexpansion(nmax,kwave,ff,1);

            % Note: slvr thinks the scatterer is at (0,0) and so
            % rad{j} thinks the origin for the expansion is (0,0)...
            % we manually set the origin to pos(j) ie the centre of the
            % scatterer
            rad{j}.origin = pos(j);
        
        end
    
        % - - - - - - - - - - - - - - - - - - - - - - - - - -
        % apply the translation addition theorem as per (35)
        % in [1]
        % - - - - - - - - - - - - - - - - - - - - - - - - - -
    
        % loop through destination scatterers
        for i=1:N
        
            % initialise expansion coefficients at the destination to zero
            reg{i} = regularzero(nmax,pos(i),kwave);
        
            % loop through all the source scatterers
            for j=1:N
                            
                if i~=j                
                    % use the translation-addition theorem to change the
                    % expansion from radiating at pos(j) to regular at
                    % pos(i)
                    reg{i} = reg{i} + regularwavefunctionexpansion(rad{j},pos(i));                
                end
                
            end
            
        end
    
        % - - - - - - - - - - - - - - - - - - - - - - - - - -
        % turn the expansions into field values on the destination
        % scatterer as per (36) in [1]
        % - - - - - - - - - - - - - - - - - - - - - - - - - -

        % loop through all the scatterer
        for j=1:N
        
            % get points on scatterer j in complex format for the incident field
            % evaluation
            qz = pos(j) + qx + 1i*qy;

            % sum the field at the Nystrom points on scatterer j
            y(:,j) = y(:,j) + 2*reg{j}.evaluate(qz);
        
        end
    
        % - - - - - - - - - - - - - - - - - - - - - - - - - -
        % finish
        % - - - - - - - - - - - - - - - - - - - - - - - - - -
        
        % turn y into a column vector
        y = y(:);

    end

    %-------------------------------------------- 
    % block Jacobi preconditioner
    %-------------------------------------------- 

    function y = preconditioner(x)

        % reshape x so that x(:,j) are the coefficients from scatterer j
        x = reshape(x,[],N);
        
        % loop through all the scatterers
        for j=1:N
            % apply the preconditioner (solve A y = x where A is the
            % Nystrom discretisation matrix, using the stored LU
            % factorisation
            y(:,j) = U \ (L \ (P*x(:,j)));
        end
        
        % turn y into a column vector
        y = y(:);

    end

end