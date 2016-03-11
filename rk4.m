function [x,y] = rk4(ystart, xstart, xend, step, derivfun)
%RK4 Implements fixed step size 4th order Runge-Kutta
%
% [x,y] = rk4(ystart, xstart, xend, step, derivfun)
%
% This function implements a 4th order Runge-Kutta integration with a fixed
% step size.  It will integrate a set of n equations.
%
% Input variables:
%
%   ystart:     n x 1 array of starting values
%
%   xstart:     scalar, starting x value
%
%   xend:       scalar, ending x value
%
%   step:       scalar, step size
%
%   derivfun:   handle to the function dy/dx = derivfun(x,y)
%
% Output variables:
%
%   x:          x values where integral was calculated
%
%   y:          y values of integrated function

% Copyright 2007 Kelly Kearney

%-----------------------------
% Setup
%-----------------------------

x = xstart:step:xend;
nx = length(x);
nstep = nx - 1;

y = zeros(length(ystart), nx);
y(:,1) = ystart;

%-----------------------------
% Loop over each time step
%-----------------------------

for istep = 1:nstep
    y(:,istep+1) = rkstep(x(istep), step, y(:,istep), derivfun);
end

%-----------------------------
% Subfunction: Compute one
% time step using Runge-Kutta
% algorithm
%-----------------------------

function y2 = rkstep(x1, h, y1, derivfun)

k1 = h * feval(derivfun, x1,       y1);
k2 = h * feval(derivfun, x1 + h/2, y1 + k1/2);
k3 = h * feval(derivfun, x1 + h/2, y1 + k2/2);
k4 = h * feval(derivfun, x1 + h,   y1 + k3);

dy = k1/6 + k2/3 + k3/3 + k4/6;
y2 = y1 + dy;


