function  err = movingMeshBurgers1D(N,tau,dt,K)
%=========================================================================%
%               MMPDE Method for 1D Burgers' Equation
%=========================================================================%
% Solve:                          
% 1D viscous Burgers' equation: 
%              ut + u ux = eps uxx (epsilon small > 0) in D = (0 1)
%              Dirichlet BCs                           on x = 0, x = 1
%
%  on an adaptive grid based on the MMPDE method. Finite differences
%  are used.
%
% Ref. W. Huang and R.D. Russell  Adaptive Moving Mesh Methods New York: 
%                                                   Springer-Verlag, 2011
% -----------------------------------------------------------------------%
% INPUTS:           - N is number of grid points (stays constant)
%                   - tau is a mesh adaptation parameter, needs to be
%                   adjusted for each problem (1 is usually fine)
%                   - dt is time step
%                   - K is number of iterations in MMPDE solve (1 or 10)
%                   
%
% OPTIONAL OUTPUTS: - xn is solution to MMPDE  
%                   - U is solution to Burgers' equation on xn
%                   - ct is computation time
%                   - err is infinity norm error to exact solution
%
% PLOTS:            - U on xn at various time steps; compare to exact solution
%                    plotted on a refined grid. 
%                   - Mesh adaptation shown on the axis
%
% CALLS:            - altSolve : MMPDE and PDE solve based on MkP
%                    procedure (see ref)
%
% Nicola M, 2017 APMA 930 CFD Final Project
%=========================================================================%
global eps 
r = cputime;
%-------------------------------------------------------------------------%
%                         Setup
%-------------------------------------------------------------------------%

eps = 1e-3;         % Coefficient of diffusive term 


endTime = 1.; numTimeSteps = endTime/dt;

% An exact solution 
Uexact = @(x,t) (0.1*exp((-x+0.5-4.95*t)/(20*eps)) + 0.5*exp((-x+0.5-0.75*t)/(4*eps)) ...
                + exp((-x+0.375)/(2*eps)))./(exp((-x+0.5-4.95*t)/(20*eps)) + ...
                exp((-x+0.5-0.75*t)/(4*eps)) + exp((-x+0.375)/(2*eps)));

%-------------------------------------------------------------------------%
%                        Initialize                                       %
%-------------------------------------------------------------------------%

% Define an initial (uniform) mesh
dx = 1/(N-1); x = (0:dx:1)';

xRef = (0:dx/1000:1)'; 

U0 = Uexact(x,zeros(N,1));



%-------------------------------------------------------------------------%
%                      Time Stepping                                      %
%-------------------------------------------------------------------------%
U1 = Uexact(x,dt);


x0 = x;
x1 = x;
for tstep = 2:numTimeSteps
    t = dt*tstep;
    
    % Generate new mesh and update solution
    [xn,U] = altSolve(x0, x1, U0,U1, tau, K, dt);
    
    % Error calc
    err(tstep) = norm(U-Uexact(xn,t),inf);
    
    % Plot snapshots
    
 if t == 10*dt
figure()
    
    plot(xn,U,'.'); 
    xlim([0,1]); ylim([-.02,1.02]); hold on;
    title(['Time = ', num2str(t)])
    plot(xRef,Uexact(xRef,t),'r');
    plot(xn,0,'b.');
    
elseif t == 0.5
figure()
    
    plot(xn,U,'.'); 
    xlim([0,1]); ylim([-.02,1.02]); hold on;
    title(['Time = ', num2str(t)])
    plot(xRef,Uexact(xRef,t),'r');
    plot(xn,0,'b.');
  
elseif t == 1
figure()
    
    plot(xn,U,'.'); 
    xlim([0,1]); ylim([-.02,1.02]); hold on;
    title(['Time = ', num2str(t)])
    plot(xRef,Uexact(xRef,t),'r');
    plot(xn,0,'b.');

end

    % Update
    U0 = U1;
    U1 = U;
    x0 = x1;
    x1 = xn;
    
    
end

err = err(end);
ct = cputime - r;

end

