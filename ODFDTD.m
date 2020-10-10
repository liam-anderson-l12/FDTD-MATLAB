function res = ODFDTD(param)
%
%   res = ODFDTD(param)
%
%   Sample 1D FDTD routine based on supplied parameter structure 
%
%   res is the results structure
%       res.Ex
%       res.Hy
%       res.Esrc
%       res.Hsrc
%       res.t
%       res.z
%

if (nargin < 1)||isempty(param)
    param = OneDParam;
end
if ~isfield(param,'fig')
    fig = 1;
else
    fig = param.fig;
end

% System parameters
z0 = param.zi(1);
zL = param.zi(end);
dmin = min(diff(param.zi));
ep0 = 8.854e-12;
mu0 = 4*pi*1e-7;
eta0 = sqrt(mu0/ep0);
c0 = 1/sqrt(ep0*mu0);
wmax = 2*pi*param.Bandwidth;
tau = 1/(2*wmax);
t0 = 6*tau;
nmax = max(sqrt(param.mur.*param.epr));
lambda_min = c0/(nmax*param.Bandwidth);

% sampling intevals
dz1= lambda_min/param.wres; 
dz2 = dmin/param.dres;
dz = min(dz1,dz2); % cell size

% grid snap to feature
N = ceil(dmin/dz); %number of cells rounded to fit smallest feature
dz = dmin/N; %cell size snapped to smallest feature
nzi = round(diff(param.zi)/dz);
nz = 1+cumsum([0 nzi]);
ns = round(param.zs/dz);

%time step 
dt=dz/(c0);

% Number of time steps
% Define Grid
nsteps = round(3*(tau/dt)+2*(zL-z0)/(dt*c0));
Nz = nz(end);
z = (0:Nz-1)*dz;
t = (0:nsteps-1)*dt;

% distributed parameters values
Ae = zeros(1,Nz);
Am = Ae;
As = Ae;
for k = (1:length(param.epr))
    Ae(nz(k):nz(k+1)-1) = param.epr(k);
    Am(nz(k):nz(k+1)-1) = param.mur(k);
    As(nz(k):nz(k+1)-1) = param.sigma(k);
end
Ae(Nz) = Ae(Nz-1);
Am(Nz) = Am(Nz-1);
As(Nz) = As(Nz-1);

% Lossy medium update coefficients
Loss = As*dt./(2*Ae*ep0);
cfe = (1 - Loss)./(1 + Loss);
cfm = 1./(1 + Loss);

% Gaussian pulse
csrc = c0/sqrt(Ae(ns)*Am(ns));
eta_e = eta0./Ae;
eta_m = eta0./Am;
Hsrc = exp(-((t - t0)/tau).^2);
Esrc = exp(-((t + dt/2 + dz/(2*csrc) - t0)/tau).^2);

%intialize fields to zero
Hx=zeros(1,Nz); Ey=zeros(1,Nz);
Hyt = zeros(nsteps,Nz); Ext = zeros(nsteps,Nz);

% run the loop
figure(fig),clf;
for nt = (1:nsteps)
    
    
    % update the magnetic field
    Hx(end) = Hx(end-1);
    Hx(1:end-1) = Hx(1:end-1) + (Ey(2:end) - Ey(1:end-1))./eta_m(1:end-1);
    Hx(ns-1) = Hx(ns-1) - Hsrc(nt)/eta_m(ns-1);
    
    % update the electric field
    Ey(1) = Ey(2);
    Ey(2:end) = cfe(2:end).*Ey(2:end) + cfm(2:end).*(Hx(2:end) - Hx(1:end-1)).*eta_e(1:end-1);
    Ey(ns) = Ey(ns) + Esrc(nt);
    
    
    plot(z,Ey,'b',z,eta0.*Hx,'r',param.zi(2)*[1 1],[-2 2],':k',param.zi(3)*[1 1],[-2 2],':k');
    ylim([-2 2]);
    title(sprintf('t = %d / %d',nt,nsteps));
   
    drawnow;
    
    % Adjust directions of fields for convenience
    Hyt(nt,:) = -Hx;
    Ext(nt,:) = Ey;
    
end
    
res.Ex = Ext;
res.Hy = Hyt;
res.Esrc = Esrc;
res.Hsrc = Hsrc;
res.t = t;
res.z = z;






