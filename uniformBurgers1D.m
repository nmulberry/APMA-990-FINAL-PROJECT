
%=========================================================================%
%               Uniform Grid for 1D Burgers' Equation
%=========================================================================%
% Implements Finite Difference scheme for the 1D Burgers' equation on a
% uniform mesh. 
% ut + u ux = eps uxx     in D = (0 1)
% Dirichlet BCs           on x=0, x=1

% SCHEME . Conservation form: ut + 0.5(u^2)x = eps uxx 

% All schemes use CNAB time stepping where diffusion term is treated implicitly.

%clear all; close all; clc;
%-------------------------------------------------------------------------%
%                           Setup 
%-------------------------------------------------------------------------%
r = cputime;

% Coefficient of diffusive term (take small and >0)
eps = 1e-3;            

% Grid
L = 1; M  = 100; dx = L/(M-1); xg = (0:dx:L)';
xr = (0:dx/500:L)'; %refined grid
dt = 1e-4; endTime = 1.; numTimeSteps = endTime/dt;

lambda = dt/dx;
sigma  = eps*dt/dx^2;

% Create system matrix, A
e = ones(M,1);
A  = spdiags([-e*sigma/2 (1+sigma)*e -e*sigma/2], -1:1, M,M);
A(1,:) = [1 zeros(1,M-1)];
A(end,:) = [zeros(1,M-1) 1];
A = sparse(A);

% Exact Solution
Uexact = @(x,t) (0.1*exp((-x+0.5-4.95*t)/(20*eps)) + 0.5*exp((-x+0.5-0.75*t)/(4*eps)) ...
                + exp((-x+0.375)/(2*eps)))./(exp((-x+0.5-4.95*t)/(20*eps)) + ...
                exp((-x+0.5-0.75*t)/(4*eps)) + exp((-x+0.375)/(2*eps)));
            
% Initial condition
Uinit = Uexact(xg,zeros(size(xg)));

%-------------------------------------------------------------------------%
%                        Solver
%-------------------------------------------------------------------------%

% Initialize
U0 = Uinit;

figure(); plot(xg,U0,'r'); hold on;
xlim([0, 1]); ylim([-1.5,1.5]);

% First time step (cheat)
U1 = Uexact(xg,dt);

% Evolve in time
for timeStep = 2:numTimeSteps
    
    t = timeStep*dt;
%------------------RHS ---------------------------------------------------% 
% 
     rhs = U1  + (1/4)*lambda*(3*([U1(end); U1(1:end-1)].^2-U1.^2) - ...
                                  ([U0(end); U0(1:end-1)].^2- U0.^2)) + ...
                 (1/2)*sigma*([U1(end); U1(1:end-1)] - 2*U1 + [U1(2:end); U1(1)]);
 

%-------------------------------------------------------------------------%  
    % Correct for BCs       
    rhs(1) = 1; rhs(end) = 0.1;
    
    % Update Solution
    U = A\rhs; 
    
    
    % Error calc 
   infErr(timeStep) = norm(U-Uexact(xg,t),inf);
    
    % Plot
if t == 10*dt
 figure()
    hold off;
    plot(xg,U,'.'); 
    xlim([0,1]); ylim([-.02,1.02]); hold on;
    title(['Time = ', num2str(t)])
    plot(xr,Uexact(xr,t),'r');
    pause(0.1);
    hold off;
elseif t == 1
figure()
    hold off;
    plot(xg,U,'.'); 
    xlim([0,1]); ylim([-.02,1.02]); hold on;
    title(['Time = ', num2str(t)])
   plot(xr,Uexact(xr,t),'r');
    pause(0.1);
    hold off;
elseif t == 2
figure()
    hold off;
    plot(xg,U,'.'); 
    xlim([0,1]); ylim([-.02,1.02]); hold on;
    title(['Time = ', num2str(t)])
    plot(xr,Uexact(xr,t),'r');
    pause(0.1);
    hold off;
end
    
    
    % Update
    U0 = U1;
    U1 = U;
    
    
    
end



ct = cputime - r;
