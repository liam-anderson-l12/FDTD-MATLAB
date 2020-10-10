function Sres = PlotODFDTD(res,param,skip)
%
%   Sres = PlotODFDTD(res,param,skip)
%
%   Plot the results of ODFDTD 
%   res is the resulty structure from ODFDTD.m
%
%   param.epr = relative permititvities of layers
%   param.mur = relative permeabilities of layers
%   param.zi = location of boundaries and interfaces (m)
%   param.zs = source location (m)
%   param.sigma = layer conductivities (S/m)
%   param.wres = wavelength spatial resolution factor 
%   param.dres = feature spatial resolution factor
%   param.Bandwidth = pulse bandwidth (one sigma in Hz)
%
%   skip = 0 (default) plot all
%   skip = 1 -> plot spectra only
%
%   Sres is a structure of spectral data
%

if (nargin < 3)||isempty(skip)
    skip = 0;
end
if (nargin < 2)||isempty(param)
    param = OneDParam();
end
if ~isfield(param,'fig')
    fig = 1;
else
    fig = param.fig;
end

%Shield, vacuum and air properties
epsilon0=8.854*10^-12; %vacuum permittivity
mu0= 4*pi*1e-7; %vacuum permeability
epsilon = epsilon0*param.epr;
mu = mu0*param.mur;
sigma = param.sigma;
eta0 = sqrt(mu0/epsilon0);

% Cell size and time stepping
t = res.t;
z = res.z;
dz = z(2)-z(1);
dt = t(2)-t(1);
[nsteps,Nz] = size(res.Ex);
nzi = round(diff(param.zi)/dz);
nz = 1+cumsum([0 nzi]);

% Initialize vectors
Ae = zeros(1,Nz);
Am = Ae;
As = Ae;
for k = (1:length(epsilon))
    Ae(nz(k):nz(k+1)-1) = epsilon(k)/epsilon0;
    Am(nz(k):nz(k+1)-1) = mu(k)/mu0;
    As(nz(k):nz(k+1)-1) = sigma(k);
end
Ae(Nz) = Ae(Nz-1);
Am(Nz) = Am(Nz-1);
As(Nz) = As(Nz-1);

% Rough calculation of terms in Poynting equation (without correcting grid locations)
S = res.Ex.*res.Hy;
U =  0.5*(epsilon0*(res.Ex.^2)*diag(Ae) + mu0*(res.Hy.^2)*diag(Am));
P = (res.Ex.^2)*diag(As);

dUdt = (U(2:end,:) - U(1:end-1,:))/(dt);
dSdx = (S(:,2:end) - S(:,1:end-1))/(dz);
C = dSdx(1:end-1,:) + dUdt(:,1:end-1) + P(1:end-1,1:end-1);

s = nextpow2(nsteps);
Ns = 4*2^s;
ExL = zeros(Ns,1);
ExL(1:nsteps) = res.Ex(:,2);
ExR = zeros(Ns,1);
ExR(1:nsteps) = res.Ex(:,end-1);
Esrc = zeros(Ns,1);
Esrc(1:nsteps) = res.Esrc;
FExL = fft(ExL);
FExR = fft(ExR);
FEsrc = fft(Esrc);
freq = (0:Ns-1)/(Ns*dt);

% Start loop
Umax = max(max(abs(dUdt)));
Smax = max(max(abs(dSdx)));
Emax = max(max(abs(res.Ex)));
Pmax = max(max(abs(P)));
Cmax = max(max(abs(C)));

if ~skip
    figure(fig);clf;
    Nplots = 5;
    for t = (2:nsteps-1)
        subplot(Nplots,1,1);
        plot(z,res.Ex(t,:),'b',z,eta0*res.Hy(t,:),'r');
        hold on;
        for k = (2:length(param.zi)-1)
            plot(z(nz(k))*[1 1],[-Emax Emax],':k');
        end
        hold off;
        xlabel('z (m)');
        ylabel('Ex (blue), \eta_0*Hy (red)'); 
        str = sprintf('t = %d/%d',t,nsteps-1);    
        title(str);
        ylim([-Emax Emax]);    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

        subplot(Nplots,1,2);
        plot(z(1:end-1),dSdx(t,:),'r',[z(1) z(end)],[0 0],'k');
        hold on;
        for k = (2:length(param.zi)-1)
            plot(z(nz(k))*[1 1],[-Smax Smax],':k');
        end
        hold off;
        xlabel('z');
        ylabel('dS/dz');    
        ylim([-Umax Umax]);      
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

        subplot(Nplots,1,3);
        plot(z,dUdt(t,:),'r',[z(1) z(end)],[0 0],'k');
        hold on;
        for k = (2:length(param.zi)-1)
            plot(z(nz(k))*[1 1],[-Umax Umax],':k');
        end
        hold off;
        xlabel('z');
        ylabel('dU/dt');    
        ylim([-Umax Umax]);   
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

        subplot(Nplots,1,4);
        plot(z,P(t,:),'r',[z(1) z(end)],[0 0],'k');
        hold on;
        for k = (2:length(param.zi)-1)
            plot(z(nz(k))*[1 1],[-Pmax Pmax],':k');
        end
        hold off;
        xlabel('z');
        ylabel('P');    
        if Pmax ~= 0
            ylim([-Pmax Pmax]);   
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

        subplot(Nplots,1,5);
        plot(z(1:end-1),C(t,:),'r',[z(1) z(end)],[0 0],'k');
        hold on;
        for k = (2:length(param.zi)-1)
            plot(z(nz(k))*[1 1],[-Cmax Cmax],':k');
        end
        hold off;
        xlabel('z');
        ylabel('C');    
        ylim([-Cmax Cmax]);   
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


        drawnow;

    end
end

indx = find(freq >= 2*param.Bandwidth);
Nmx = indx(1);
figure(fig+1),clf;
subplot(3,1,1);
R = abs(FExL(1:Nmx)./FEsrc(1:Nmx)).^2;
semilogy(freq(1:Nmx)*1e-9,R,'b');
xlim([freq(1) freq(Nmx)]*1e-9);
xlabel('Frequency (GHz)');
ylabel('R');
subplot(3,1,2);
T = abs(FExR(1:Nmx)./FEsrc(1:Nmx)).^2;
semilogy(freq(1:Nmx)*1e-9,T,'b');
xlim([freq(1) freq(Nmx)]*1e-9);
xlabel('Frequency (GHz)');
ylabel('T');
subplot(3,1,3);
A = 1 - (R + T);
plot(freq(1:Nmx)*1e-9,A,'b',[freq(1) freq(Nmx)]*1e-9,[0 0],':k');
axis([freq(1)*1e-9 freq(Nmx)*1e-9 -0.1 1]);
xlabel('Frequency (GHz)');
ylabel('A = 1-(R+T)');


Sres.FExR = FExR;
Sres.FExL = FExL;
Sres.FEsrc = FEsrc;
Sres.freq = freq;


