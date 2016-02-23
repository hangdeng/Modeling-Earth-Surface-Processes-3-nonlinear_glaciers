% 1d glacier model
% hang deng
% used RSA's written code
% erosion rate too small, therefore, the erorate was multipled by 500(?)
% better erosion equation can be found from Macgregor et al, geology paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
figure(1)
clf


%% initialize
% parameter
% Ice rheology and material properties
rho1 = 917; % ice density [kg m^-3]
A = 2.1e-16; % Pa-3 yr-1

% Topography
zbMax = 4000; % maximum bedrock elevation [m]
zbSlope = -0.03; % bedrock surface slope [m/m]
w0 = 1000; % valley width

% Mass balance / climate
zEla = 3700; % elevation of the equilibirum line elevation [m]
E = 400; % deltaE
dbdz = 0.01; % mass balance lapse rate [m y^-1 m^-1] (usually 0.01 m/y/m or so)
bCap = 1; % maximum mass balance [m y^-1]


% Model domain
% Space
dx = 100; % spatial step [m]
xMax = 30000; % end of spatial domain [m]

% Time
dt = 0.0025; % timestep [y] (this needs to be small [2.5e-3 works] for numerical stability)
tmax = 750; % end of temporal domain [y]


% constants
g = 9.81; % acceleration due to gravity [m s^-2]

% plotting

nplots = 50;
tplot = tmax/nplots;

%% loop and run

% Initializing space and time vectors
x = dx/2:dx:xMax-(dx/2); % x-coordinates of node centers [m]
xedge = 0:dx:xMax; % x-coordinate of node edges [m]
t = 0:dt:tmax;

% harmonically ELA
%zEla = Ela + E .* sin(0.01*t);

zb = zbMax + zbSlope*x; % bedrock elevation [m]
H = zeros(size(x)); % ice thickness
erosion = zeros(size(x));
zs = zb+H-erosion; % ice surface elevation [m]
zbed = zb -erosion;
b = dbdz*(zs-zEla); % mass balance profile [m y^-1]

for i = 1:length(t)
    % Recalculate mass balance
    b = dbdz*(zs-zEla); % local net balance calculated at cell centers at ice surface
    b = min(b,bCap); % limit to prescribed maximum mass balance

    Hedge = H(1:end-1)+0.5*diff(H); % interpolates ice thickness to cell edges
    alpha = diff(zs)/dx; % slope of ice surface calculated at cell edges

    Q = (2*A/5).*((rho1*g*sin(-alpha)).^3).*(Hedge.^5); % ice discharge due to internal deformation
    Q =[0 Q 0]; %takes care of boundary conditions - no flux in at top and no flux out at bottom
    % other option: Q(find(Q<0))=0;
    
    % erosion module added
    erorate = 500*(diff(Q)/dx)/w0; % get into erosion, equal to coefficient*flux/valley width
    erosion = erosion + erorate.*dt; 
    dHdt = b - (diff(Q)/dx); % time rate of change of ice thickness [m y^-1]
    H = H + (dHdt*dt); % update ice thickness
    H = max(H,0); % can't have negative thickness
    
    erosion(find(erosion<0))=0; % erorate > 0, make sure erosion > 0
    % then use '-' to simulate erosion
    zs = zb+H-erosion; %updates surface elevation
    zbed = zb-erosion;


    % now for some plotting
    if rem(t(i),tplot)==0
        disp(['Time: ' num2str(t(i))])
        figure(1)
        plot(x/1000,zbed,'m-','linewidth',2)
        axis([0 x(end)/1000 3000 5000])
        hold on
        plot(x/1000,zs,'c-','linewidth',1)
        xlabel('Down-glacier distance [km]','fontsize',18)
        ylabel('Elevation[m]','fontsize',18)
        legend('Bedrock after Erosion','Ice surface','location','northwest')
        set(gca,'fontsize',14)
        pause(0.01)

    end

end
