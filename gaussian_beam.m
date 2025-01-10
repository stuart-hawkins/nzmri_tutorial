% Gaussian beam object.
%
%  p = gaussian_beam(theta,x0,w0,k) returns a Gaussian beam object p with
%  wavenumber k, waist located at x0, waist radius w0 and direction
%  exp(1i*theta).
%
% Also:
%
%   f = p.evaluate(z) returns the values f of the beam at points z.
%
%   f = u.evaluate(z,mask) returns the values f of the beam at
%   points z for which mask==1 and NaN elsewhere.
%
%   [dx,dy] = u.evaluateGradient(z) returns dx and dy the partial
%   derivatives of the beam in the x and y directions respectively
%   at the points z.
%
%   [dx,dy] = u.evaluateGradient(z,mask) returns dx and dy the partial
%   derivatives of the beam in the x and y directions respectively
%   at the points z for which mask==1 and NaN elsewhere.
%
% Note: in the functions above vectors in the plane are represented by
% complex numbers.
%
% References:
%
%  [1] https://doc.comsol.com/6.1/doc/com.comsol.help.roptics/roptics_ug_optics.6.62.html
%
% See also: point_source, incident, plane_wave.
%
% Stuart C. Hawkins - 10 January 2025

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


classdef gaussian_beam < incident

    properties
        direction
        kwave
        origin
        w0
        zR
    end

    methods

        %-----------------------------------------------
        % constructor
        %-----------------------------------------------

        function self = gaussian_beam(varargin)

            if nargin==1

                % copy from an existing gaussian_beam....

                if ~isa(varargin{1},'gaussian_beam')

                    error('Single argument must be a gaussian_beam')

                end

                % set wavenumber
                self.kwave = varargin{1}.kwave;

                % set incident direction
                self.direction = varargin{1}.direction;

                % set origin
                self.origin = varargin{1}.origin;

                % set waist radius
                self.w0 = varargin{1}.w0;

                % set Rayleigh range
                self.xR = varargin{1}.xR;

            else

                % set wavenumber
                self.kwave = varargin{4};

                % set incident direction
                self.direction = exp(1i*varargin{1});

                % set origin
                self.origin = varargin{2};

                % set waist radius
                self.w0 = varargin{3};

                % compute wavelength
                lambda = 2*pi/self.kwave;

                % set Rayleigh-range
                self.zR = pi*self.w0^2/lambda;

            end

        end

        %-----------------------------------------------
        % return vector of coefficients for scatterer at
        % centre
        %-----------------------------------------------

        function cof = get_coefficients(self,centre,nmax)

            error('This operation is not yet supported for gaussian_beam')

        end

        %-----------------------------------------------
        % evaluate method
        %-----------------------------------------------

        function val = evaluate(self,points,mask)

            % intialize return array
            val=zeros(size(points));

            % apply mask if necessary
            if nargin>2
                points=points(mask);
            end
           
            % compute (u,v) which are the projections of points into the
            % beam direction and its normal
            u = real((points-self.origin)*conj(self.direction));
            v = real((points-self.origin)*conj(1i*self.direction));

            % evaluate the beam value according to [1]
            % Note: we have conjugated the expression in [1] so that the
            % beam travels in the positive u-direction for exp(-i omega t)
            % time dependence
            w = self.w(u);
            f = sqrt(self.w0./w).*exp(-v.^2./w.^2+1i*self.kwave*u ...
                +0.5i*self.kwave*v.^2./self.R(u) - 1i*self.eta(u) );

            % insert values into the return array
            if nargin>2
                val(mask)=f;
            else
                val=f;
            end

        end

        %-----------------------------------------------
        % evaluate gradient
        %
        % This is useful for eg Neumann BCs. Note that we
        % cannot use complex numbers to represent vectors
        % in this case because the components of the vector
        % are, in general, complex.
        %-----------------------------------------------

        function [dx,dy] = evaluateGradient(self,points,mask)

            % initialize return values
            dx = zeros(size(points));
            dy = zeros(size(points));

            % apply mask if necessary
            if nargin>2
                points=points(mask);
            end

            % compute (u,v) which are the projections of points into the
            % beam direction and its normal
            u = real((points-self.origin)*conj(self.direction));
            v = real((points-self.origin)*conj(1i*self.direction));

            % evaluate subfunctions
            w = self.w(u);
            dw = self.dw(u);
            R = self.R(u);
            dR = self.dR(u);

            % evaluate exponential... will use this in the partial
            % derivatives
            ff = exp( -v.^2./w.^2+1i*self.kwave*u ...
                +0.5i*self.kwave*v.^2./R - 1i*self.eta(u) );

            % compute partial derivatives with respect to u and v
            dfdu = -0.5*sqrt(self.w0)./w.^(3/2).*dw.*ff + ...
                sqrt(self.w0./w).*( 2*v.^2./w.^3.*dw + 1i*self.kwave ...
                -0.5i*self.kwave*v.^2./R.^2.*dR - 1i*self.deta(u) ).*ff;

            dfdv = sqrt(self.w0./w).*( -2*v./w.^2 + 1i*self.kwave*v./R ).*ff;           
            
            % compute partial derivatives of the projection onto the beam
            % direction and its normal
            dudx = real(self.direction);
            dudy = imag(self.direction);
            dvdx = -imag(self.direction);
            dvdy = real(self.direction);

            % use the chain rule to get the derivatives with respect to x
            % and y
            dfdx = dfdu*dudx + dfdv*dvdx;
            dfdy = dfdu*dudy + dfdv*dvdy;

            % insert values into the return array
            if nargin>2
                fx(mask)=dfdx;
                fy(mask)=dfdy;
            else
                fx=dfdx;
                fy=dfdy;
            end

            % insert values into the return array
            if nargin>2
                dx(mask) = fx;
                dy(mask) = fy;
            else
                dx = fx;
                dy = fy;
            end

        end

        %-----------------------------------------------
        % various subfunctions used to evaluate the beam
        %-----------------------------------------------

        % Note: prefix d indicates derivate eg dw is derivative of w

        function val = w(self,x)

            val = self.w0 * sqrt(1+(x/self.zR).^2);

        end

        function val = dw(self,x)

            val = self.w0 * (x/self.zR.^2)./sqrt(1+(x/self.zR).^2);

        end

        function val = R(self,x)

            val = x .* (1+(self.zR./x).^2);

        end

        function val = dR(self,x)

            val =  (1+(self.zR./x).^2) - 2*self.zR^2./x.^2;

        end

        function val = eta(self,x)

            val = 0.5*atan(x/self.zR);

        end

        function val = deta(self,x)

            val = 0.5/self.zR*1./(1+(x/self.zR).^2);

        end

    end % end methods

end