function [x,f] = altSolve(x0, x1, f0, f1, tau, K , dt)
% Get mesh based on MkP procedure; solve Burgers' equation on this mesh.
%       Spatial Discretization
%              - MMPDE is discretized on computational grid using a
%                finite difference scheme
%              - PDE is discretized on mesh  using a finite difference
%                scheme
%       Time Stepping:
%              - Backwards Euler for MMPDE
%              - CNAB for physical PDE
% INPUTS: - f0 is solution at time step n-1, defined on mesh x0
%         - f1 is solution at time step n, defined on mesh x1
%         - tau  is MMPDE paramater
%         - K is number of iterations in MMPDE solve
%         - dt is time step
%
% OUTPUTS: - x is new mesh
%          - f is new solution 
% ------------------------------------------------------------------------%

global eps
N      = length(x0);
dtn    = dt/K;

dx      = 1/(N-1);

%-------------------MESH DENSITY FN---------------------------------------%

D       = D2U(x1,f1,N);
alpha   = max(1, (1/8)*sum(diff(x1).*(abs(D(2:N)).^(2/3) + abs(D(1:N-1)).^(2/3)))^3);
rho     = (1+ (1/alpha)*abs(D).^2).^(1/3); 


rhon = rho;
xn = x1;

for k = 1:K

    % Smooth mesh density and interpolate to new grid
    rho = interp1(x1,rhon,xn,'linear');
    for jj = 1:4
        rho     = [(1/2)*(rho(1) + rho(2)); ...
                   (1/4)*rho(1:N-2)+(1/2)*rho(2:N-1)+(1/4)*rho(3:N); ...
                   (1/2)*(rho(N-1)+rho(N))];
    end

    %---------------------SOLVE MMPDE-----------------------------------------%

A = zeros(N,N);
A(1,1) = 1;
    for ii = 2:N-1
        rhof = (rho(ii+1) + rho(ii))/2;
        rhob = (rho(ii-1) + rho(ii))/2;

        A(ii,ii+1) = -dtn*rhof/(dx^2*tau*rho(ii));
        A(ii,ii-1) = -dtn*rhob/(dx^2*tau*rho(ii));
        A(ii,ii) = 1-A(ii,ii+1)-A(ii,ii-1);

    end
A(N,N) = 1;

x = A\xn;    x = sort(x);     

xn = x;

end

%---------------------SOLVE BURGERS---------------------------------------%

% Build system matrix (needs to be updated at every iteration)
s = (x - x1)./dt; 
A = zeros(N,N);
rhs = zeros(size(f0));

rhs(1) = f1(1);
A(1,1) = 1;

for ii = 2:N-1

% Build system matrix using grid at time n+1
    dx0 = 1/(x(ii+1) - x(ii-1));
    dxb = 1/(x(ii) - x(ii-1));
    A(ii,ii-1) = -eps*dt*dxb*dx0;
    dxf = 1/(x(ii+1) - x(ii));
    A(ii,ii+1) = -eps*dt*dxf*dx0; 
    A(ii,ii) = 1 - A(ii,ii-1) - A(ii,ii+1);

% Build RHS using grid at time n
    dx0 = 1/(x1(ii+1) - x1(ii-1));
    dxb = 1/(x1(ii) - x1(ii-1));
    dxf = 1/(x1(ii+1) - x1(ii));

    rhs(ii) = f1(ii) + (3/4)*dt*dxb*(f1(ii-1)^2 - f1(ii)^2) + ...
                    (3/2)*dt*dx0*(s(ii))*(f1(ii+1) - f1(ii-1)) + ...
              (1/2)*eps*dt*dx0*(dxf*(f1(ii+1) - f1(ii)) - dxb*(f1(ii) - f1(ii-1)));

% Build RHS using grid at time n-1    
    dxb = 1/(x0(ii) - x0(ii-1));
    dx0 = 1/(x0(ii+1) - x0(ii-1));
    dxf = 1/(x0(ii+1) - x0(ii));

    rhs(ii) = rhs(ii) - (1/4)*dt*dxb*(f0(ii-1)^2 - f0(ii)^2) - ...
                        (1/2)*dt*dx0*(s(ii))*(f0(ii+1) - f0(ii-1));

end

A(N,N) = 1;
rhs(N) = f1(N);

A = sparse(A);
f = A\rhs;

end

function uxx = D2U(x,u,N)

% Interior Points 2:N-1

xFor = x(3:N)-x(2:N-1);
xBack = x(2:N-1)-x(1:N-2);
xCent = x(3:N)-x(1:N-2);

uFor = u(3:N)-u(2:N-1);
uBack = u(2:N-1)-u(1:N-2);

uxx = (2./xCent).*(uFor./xFor - uBack./xBack);

% Boundary Points 1,N

uxx1 = 2*((x(2)-x(1))*(u(3) - u(1)) - (x(3)-x(1))*(u(2)-u(1))) / ...
          ( (x(3)-x(1))*(x(2)-x(1))*(x(3)-x(2)) );
      
uxxN = 2*((x(N-1)-x(N))*(u(N-2) - u(N)) - (x(N-2)-x(N))*(u(N-1)-u(N))) / ...
          ( (x(N-2)-x(N))*(x(N-1)-x(N))*(x(N-2)-x(N-1)) );
      
% Output
uxx = [uxx1;uxx;uxxN];

end


