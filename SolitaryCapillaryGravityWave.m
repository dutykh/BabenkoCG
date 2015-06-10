function [zs, ws, fs, SWP, W, F, P, A] = SolitaryCapillaryGravityWave (Fr, Bo, eta0, Z, PF)
% SolitaryCapillaryGravityWave: This function computes 
% the steady irrotational surface solitary (classical 
% and generalized) capillary-gravity wave solutions of 
% the full Euler equations (homogeneous, incompressible 
% and perfect fluids). The wave is defined by its initial 
% Froude and Bond numbers (Fr, Bo) and the result is about 
% twelve digits accurate. The method works for all but the 
% highest waves.
% 
% NOTA BENE: The execution of this function requires the 
% presence of the Optimization Toolbox, namely we use the 
% fsolve() function of this toolbox to solve the nonlinear
% nonlocal Babenko equation.
%
% SYNOPSIS:
% [zs,ws,fs,SWP] = SolitaryCapillaryGravityWave(Fr, Bo);
% [zs,ws,fs,SWP,W,F,P,A] = SolitaryCapillaryGravityWave(Fr, Bo, eta);
% [zs,ws,fs,SWP,W,F,P,A] = SolitaryCapillaryGravityWave(Fr, Bo, eta, Z);
% [zs,ws,fs,SWP,W,F,P,A] = SolitaryCapillaryGravityWave(Fr, Bo, eta, Z, 1);
% 
% INPUT:
% Fr : Froude number (i.e., c/sqrt(gd), c.f. Note 1 below).
% Bo : Bond number (i.e. \tau/(g*d^2), c.f. Note 1 below).
% eta: Initial guess of the solution (vector of N real numbers)
% Z  : Complex abscissa where fields are desired inside the fluid (default Z = []).
%      Z should be strictly below the surface, i.e., -1 <= imag(Z) < eta(real(Z))
%      y = eta(x) being the equation of the free surface.
% 
% Have a look at the m-file for more details.
% 
% OUTPUT (dimensionless quantities):
% zs   : Complex abscissa at the surface, i.e., x + i*eta.
% ws   : Complex velocity at the surface, i.e., u - i*v.
% fs   : Complex potential at the surface, i.e., phi + i*psi.
% SWP  : Solitary Wave Parameters, i.e.
%        SWP(1) = wave height, max(eta) - min(eta)
%        SWP(2) = wave mass
%        SWP(3) = circulation
%        SWP(4) = impulse
%        SWP(5) = kinetic energy
%        SWP(6) = potential gravity energy
%        SWP(7) = potential capillary energy
% W    : Complex velocity in the bulk at abscissas Z.
% F    : Complex potential in the bulk at abscissas Z.
% P    : Pressure in the bulk at abscissas Z.
% A    : Complex acceleration in the bulk at abscissas Z (A = dW / dt).
% 
% NOTES:
% 1- Computations are performed in dimensionless units such that rho = g = d = 1,
%    where rho is the fluid density, g is the acceleration due to gravity and d 
%    is the constant water depth. It means that all lengths are scaled by d,
%    accelerations by g, speeds by sqrt(g*d), times by sqrt(d/g) and stresses 
%    by rho*g*d.
% 2- The numerical scheme is based on the Levenberg-Marquardt algorithm to solve
%    a system of nonlinear equations resulting from the Babenko equation
%    discretization using a pseudo-spectral collocation method.
% 3- The solution is obtained in parametric form resulting from a conformal 
%    mapping into a strip.
% 4- The algorithm is efficient for all but the highest waves.
% 5- W, F, P and A are column vectors corresponding to the shape of Z(). The
%    velocities are computed in the frame of reference where the fluid is
%    at rest in the far field.
% 6- For Z very close to the free surface, the fields may be inaccurate.
% 
% REMARKS:
% If the method does not converge for some reasonable values 
% Fr, Bo parameters, there are several possible solutions:
%		1. Change the initial guess
%       2. Change the length of the domain
% 		3. Employ the natural continuation in one of the
%		parameters from a fully converged solution
% The numerical computation of capillary-gravity waves is often
% a tricky matter and some numerical experience is needed to obtain 
% new and non-trivial solutions. For some examples, please consult
% the REFERENCE given below.
% 
% Author 1: Denys Dutykh, CNRS & LAMA UMR 5127
%			University Savoie Mont Blanc
% E-mail  : Denys.Dutykh@univ-savoie.fr
% URL     : http://www.denys-dutykh.com/
%
% Author 2: Didier Clamond, LJAD
% 			University of Nice - Sophia Antipolis, France
% E-mail  : Didier.Clamond@gmail.com
% 
% Author 3: Angel Duran, University of Valladolid, Spain
% E-mail  : Angel.Duranmartin0@gmail.com
% 
% REFERENCE: D. Clamond, D. Dutykh & A. Duran. A plethora 
% of generalised solitary gravity-capillary water waves.
% Submitted, 2015
% 	  https://hal.archives-ouvertes.fr/hal-01081798/

if (nargin < 4)
	Z = [];
	PF = 0;
end % if ()

if (nargin < 5)
	PF = 0;
end % if ()

%%% Physical parameters:
g   = 1.0;          % gravity acceleration
d   = 1.0;          % undisturbed water depth
c0  = sqrt(g*d);    % linear gravity wave speed
c   = Fr/c0;        % solitary wave speed
c2  = c^2;          % Fr^2
sig = Bo/(g*d*d);   % unscaled surface tension

if (c > 1.29421)
   error('The Froude number must be in the interval Fr <= 1.29421');
end % if ()

if (sig < 0.0)
   error('The Bond number must be positive');
end % if ()

if ((c > 0) && (sig > 1.0/3.0))
	error('If Fr > 0 than Bo must be < 1/3');
end % if ()

if ((c < 0) && (sig < 1.0/3.0))
	error('If Fr < 0 than Bo must be > 1/3');
end % if ()

% Decay of the gravity solitary wave
kd = fzero(@(kd) sinc(kd)-c2*cos(kd), [0; 1.5]); % trend parameter
% To be eventually modified:
L   = 17.0*log(10)/kd;	% half-length of the computational domain
N   = 1024;     		    % number of Fourier modes (must be even)
tol = 1e-12;    		    % iterations stopping criterium
dxi = 2*L/N;            % distance between two points in conformal space
xi  = (1-N/2:N/2)'*dxi; % conformal space discretization
k   = [0:N/2-1, -N/2:-1]'*pi/L; % wavenumbers

%%% Initial guess for the solution (classical Serre equations):
if ((nargin < 3) || (isempty(eta0)))
    eta0 = (c2 - 1.0)*sech(0.5*kd*xi).^2;
end % if()

%%% Definition of some pseudo-differential operators:
T   = 1i*tanh(k*d);
C   = k.*coth(k*d); C(1) = 1/d;     % nonlocal C-operator in the Babenko equation
Ci  = tanh(k*d)./(k*d); Ci(1) = 1;  % inverse C operator
Cnl = -1i*coth(k*d); Cnl(1) = 0;	% nonlocal int-C-operator

%%% The linear operator in Babenko equation:
L = c2 - Ci - sig*k.*tanh(k*d);

%%% Finally, we run the nonlinear solver:
opts = optimset('Display', 'iter', 'FinDiffType', 'forward',...
    'MaxFunEvals', 200*N, 'Algorithm', {'levenberg-marquardt', 0.05},...
    'TolX', tol, 'TolFun', tol);
[eta, fval, exitflag, output] = fsolve(@Babenko, eta0, opts);

switch exitflag
	case 1
		fprintf('Function converged to a solution.\n');
	case 2
		fprintf('Change in x was smaller than the specified tolerance.\n');
	case 3
		fprintf('Change in the residual was smaller than the specified tolerance.\n');
	case 4
		fprintf('Magnitude of search direction was smaller than the specified tolerance.\n');
	case 5
		fprintf('Number of iterations exceeded MaxIter or number of function evaluations exceeded MaxFunEvals.\n')
	case -1
		fprintf('Output function terminated the algorithm.\n');
	case -2
		fprintf('Algorithm appears to be converging to a point that is not a root.\n');
	case -3
		fprintf('Trust region radius became too small.\n');
	case -4
		fprintf('Line search cannot sufficiently decrease the residual along the current search direction.\n');
end % switch

% Post processing.
Ceta   = real(ifft(C.*fft(eta)));                 % C(eta)
dexi   = real(ifft(1i*k.*fft(eta)));              % d eta / d xi
SWP(1) = max(eta) - min(eta);                     % Wave height
Ceta1  = 1.0 + Ceta;
SWP(2) = dxi*eta'*Ceta1;                          % mass
SWP(3) = c*dxi*sum(Ceta);                         % circulation
SWP(4) = c*SWP(2);                                % impulse
SWP(5) = 0.5*c2*dxi*eta'*Ceta;                    % kinetic energy
SWP(6) = 0.5*dxi*(eta.^2)'*Ceta1;                 % potential gravity energy
% potential capillary energy:
SWP(7) = sig*dxi*sum(sqrt(Ceta1.^2 + dexi.^2) - Ceta1);

% Physical variables at the surface.
etaMean = mean(eta);
xs   = (1 + etaMean)*xi + real(ifft(Cnl.*fft(eta - etaMean))); 
phis = c*(xs - xi);
qs   = (1+Ceta).^2 + dexi.^2;
us   = c - c*(1+Ceta)./qs;
vs   = -c*dexi./qs;

% Output complex variables at the surface.
zs = xs + 1i*eta;
ws = us - 1i*vs;
fs = phis + 1i*eta*c;

fprintf('+-------------------------------------------------+\n');
fprintf('| Froude Number           = %15.14f      |\n', c);
fprintf('| Bond Number             = %15.14f      |\n', sig);
fprintf('| Wave height             = %15.14f      |\n', SWP(1));
fprintf('| Mass                    = %15.14f      |\n', SWP(2));
fprintf('| Circulation             = %15.14f      |\n', SWP(3));
fprintf('| Impulse                 = %15.14f      |\n', SWP(4));
fprintf('| Kinetic Energy          = %15.14f      |\n', SWP(5));
fprintf('| Potential Gr. Energy    = %15.14f      |\n', SWP(6));
fprintf('| Potential Cap Energy    = %15.14f      |\n', SWP(7));
fprintf('|                                                 |\n');
fprintf('| Convergence achieved in %05.0f iterations.       |\n', output.iterations);
fprintf('| Error between two latest iterations: %5.4e |\n', output.stepsize);
fprintf('| Residual                           : %5.4e |\n', norm(fval, inf));
fprintf('+-------------------------------------------------+\n');

% Output at desired locations in the bulk.
% The code below is not fully vectorize in order to avoid building a huge 
% matrix when Z is large.
if isempty(Z)==0,
   Z = Z(:);
   LZ = length(Z);
   W = zeros(LZ,1); F=W; P=W; dW=W;
   dzsm1 = Ceta + 1i*dexi;
   for n=1:LZ,
      W(n)  = sum(dzsm1./(zs-Z(n)) - conj(dzsm1)./(conj(zs)-2i-Z(n))); 
      W(n)  = 1i*0.5*c/pi*W(n)*dxi;
      dW(n) = sum(dzsm1./(zs-Z(n)).^2 - conj(dzsm1)./(conj(zs)-2i-Z(n)).^2); 
      dW(n) = 1i*0.5*c/pi*dW(n)*dxi;
      F(n) = sum(dzsm1.*log((zs+1i)./(zs-Z(n))) - conj(dzsm1.*log((zs+1i)./(zs+2i-conj(Z(n)))))); 
      F(n) = 1i*0.5*c/pi*F(n)*dxi;
   end
   P = c*real(W) - 0.5*abs(W).^2 - imag(Z); % pressure
   A = dW.*(conj(W) - c);                   % acceleration
end

if (PF)
   Lc = 2*SWP(2)/SWP(1); % characteristic length
   
   subplot(2,2,1)
   plot(xs, eta, 'b-','LineWidth',1)
   title('Free surface elevation', 'interpreter', 'latex', 'fontsize', 12)
   xlabel('$x\ /\ d$', 'interpreter', 'latex', 'fontsize', 14)
   ylabel('$\eta\ /\ d$', 'interpreter', 'latex', 'fontsize', 14)
   axis([-Lc Lc 1.05*min(eta) 1.05*SWP(1)])

   subplot(2,2,2)
   plot(xs, phis, 'b-','LineWidth',1)
   title('Velocity potential', 'interpreter', 'latex', 'fontsize', 12)
   xlabel('$x\ /\ d$', 'interpreter', 'latex', 'fontsize', 14)
   ylabel('$\phi\ /\ d\ \sqrt{\, g\/d\ }$', 'interpreter', 'latex', 'fontsize', 14)
   axis([-Lc Lc 1.05*min(phis) 1.05*max(phis)])
   
   subplot(2,2,3)
   plot(xs,us, 'b-','LineWidth',1)
   title('Horizontal velocity', 'interpreter', 'latex', 'fontsize', 12)
   xlabel('$x\ /\ d$', 'interpreter', 'latex', 'fontsize', 14)
   ylabel('$u\ /\ \sqrt{\, g\/d\ }$', 'interpreter', 'latex', 'fontsize', 14)
   axis([-Lc Lc 1.05*min(us) 1.05*max(us)])
   
   subplot(2,2,4)
   plot(xs,vs, 'b-','LineWidth',1)
   title('Vertical velocity', 'interpreter', 'latex', 'fontsize', 12)
   xlabel('$x\ /\ d$', 'interpreter', 'latex', 'fontsize', 14)
   ylabel('$v\ /\ \sqrt{\, g\/d\ }$', 'interpreter', 'latex', 'fontsize', 14)
   axis([-Lc Lc 1.05*min(vs) 1.05*max(vs)]);

   drawnow;
end % if (PF)

function Eq = Babenko(eta)
    eta2     = eta.*eta;
    eta2_hat = fft(eta2);
    eta_hat  = fft(eta);    
    eta_xi   = real(ifft(1i*k.*eta_hat));
    Ceta     = real(ifft(C.*eta_hat));
    Ceta1    = 1.0 + Ceta;
    root     = sqrt(Ceta1.*Ceta1 + eta_xi.*eta_xi);

    Eq = real(ifft(C.*(c2*eta_hat - 0.5*g*eta2_hat + sig*fft(1 - Ceta1./root))));
    Eq = Eq - g*eta.*Ceta1 + sig*real(ifft(1i*k.*fft(eta_xi./root)));
end % Babenko ()

% SINC(X) is SIN(X)/X  with the singularity at zero removed.
function y = sinc(x)
	% Drea Thomas (1992)
  	z     = (x ~= 0);
  	y     = x;
  	y(z)  = sin(x(z))./x(z);
  	y(~z) = ones(sum(sum(~z)),1);
end % sinc()

end % SolitaryCapillaryGravityWave ()