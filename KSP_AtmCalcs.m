function [t,x] = KSP_AtmCalcs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculate optimal atmospheric ascend profile and/or maximum altitude and
% speed obtainable with a winged spaceplane with KSP's (0.24.2) Ion Engine.
%
% Five forces are accounted for: lift, drag, thrust, gravity, centrifugal.
%
% 1) Lift:        FL = [-vz;vy]*L*Cl*P(z)
% 2) Drag:        FD = -|v|*[vy;vz]*0.0049*Cd*m*P(z)
% 3) Thrust:      FT = T/|v|*[vy;vz]
% 4) Gravity:     FG = [0;-1/(R+z)^2]*G*m*M
% 5) Centrifugal: FC = [0;vy^2]*m/(R+z)
%
% ODE45 is used to solve the 2 coupled (nonlinear) second-order differential
% equations, by first re-writing them into a system of 4 coupled (nonlinear)
% first-order differential equations.
%
% Author: Kay Steinkamp, Dec 2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Part Specs as of 0.24.2 
% -------------------------------------------------------
IonEngine.thrust = 2;
IonEngine.mass = 0.25;
IonEngine.TWR = IonEngine.thrust/IonEngine.mass/9.81;
IonEngine.XpS = 0.485;
IonEngine.EpS = 8.729;
IonEngine.Isp = 4200;

XenonS.wetmass = 0.07;
XenonS.amount = 400;
XenonS.drymass = 0.03;
XenonL.wetmass = 0.12;
XenonL.amount = 700;
XenonL.drymass = 0.05;

Swept.mass = 0.05;
Swept.maxdrag = 0.6;
Swept.lift = 1.6;

DeltaD.mass = 0.02;
DeltaD.maxdrag = 0.6;
DeltaD.lift = 0.7;

CtrlS.mass = 0.01;
CtrlS.maxdrag = 0.5;
CtrlS.lift = 0.5;

Okto2.mass = 0.04;


% Planet characteristics as of 0.24.2 
% -------------------------------------------------------
Kerbin.R = 6e5;       % m
Kerbin.M = 5.2916e22; % kg
Kerbin.P0 = 1;        % atm
Kerbin.H = 5000;      % m
Kerbin.Vr = 2*pi*Kerbin.R/2.16e4; % Rotational velocity at surface

Eve.R = 7e5; 
Eve.M = 1.2244e23;
Eve.P0 = 5; 
Eve.H = 7000; 
Eve.Vr = 2*pi*Eve.R/8.05e4; 

Laythe.R = 5e5; 
Laythe.M = 2.9398e22;
Laythe.P0 = 0.8; 
Laythe.H = 4000; 
Laythe.Vr = 2*pi*Laythe.R/5.298e4;

Duna.R = 3.2e5; 
Duna.M = 4.5155e21;
Duna.P0 = 0.2; 
Duna.H = 3000; 
Duna.Vr = 2*pi*Duna.R/6.552e4;


% Fixed parameters
% -------------------------------------------------------
planet = 'Duna';

global R M G P0 H Vr minDrag restDrag
R = eval([planet '.R;']);
M = eval([planet '.M;']);
P0 = eval([planet '.P0;']);
H = eval([planet '.H;']);
Vr = eval([planet '.Vr;']);

G = 6.673e-11; % grav. const.
minDrag = 0.02; % minimum drag coefficient for wing parts
restDrag = 0.2;  % drag of rest of plane (non-wing parts)
wmaxDrag = Swept.maxdrag; % max drag of wing-like parts
cmaxDrag = CtrlS.maxdrag; % max drag of controlsurface-like parts


% Ship design
% (ignore massless parts)
% -------------------------------------------------------
global Cd Cl_w Cl_s L_w L_s wingDrag ctrlDrag m0 mf mwings mctrls tburn T AoIs

nIE = 2;  % number of ion engines
nXS = 1;  % number of small Xenon tanks
nXL = 4;  % number of large Xenon tanks
nSW = 2;  % number of swept wings;
nDD = 0;  % number of Delta Deluxe winglets
nCS = 5;  % number of small control surfaces
AoIs = 2; % Angle of Incidence for the control surfaces

mwings = 1e3*(nSW*Swept.mass + nDD*DeltaD.mass); % mass of all wings (maxDrag=0.6, wing-like L)
mctrls = 1e3*(nCS*CtrlS.mass); % mass of all control surfaces (maxDrag=0.5, ctrlS-like L)
m0 = mwings + mctrls + 1e3*(nIE*IonEngine.mass + nXS*XenonS.wetmass ...
                            + nXL*XenonL.wetmass + Okto2.mass); % kg
mf = m0 - 1e3*(nXS*(XenonS.wetmass-XenonS.drymass) ...
              + nXL*(XenonL.wetmass-XenonL.drymass)); %kg
tburn = (nXS*XenonS.amount + nXL*XenonL.amount)/(nIE*IonEngine.XpS); % s
T = nIE*IonEngine.thrust; % kN
Cl_w = nSW*Swept.lift + nDD*DeltaD.lift; % wings lift coeff. 
Cl_s = nCS*CtrlS.lift;                   % control surfaces lift coeff. 
wingDrag = mwings*wmaxDrag; % AoA-dependent drag coeff. for wings
ctrlDrag = mctrls*cmaxDrag; % AoA-dependent drag coeff. for ctrl surfaces
Cd = []; % total (mass-averaged) drag coeff., computed in subroutine
L_w = [];  % AoA-dependent lift function for wings, computed in subroutine
L_s = [];  % AoA-dependent lift function for ctrls, computed in subroutine


% Initial values
% -------------------------------------------------------
global AoW AoA AoV AoAbounds tprev
t0 = 0;               % initial time
tf = 1.1*tburn;       % final time
tprev = 0;            % will contain previous ODE timestep
x10 = 0;              % initial horizontal position y (Note, x1=y)
x20 = 12;             % initial horizontal velocity yd (Note, x2=yd)
x30 = 70;             % initial vertical position z (Note, x3=z)
x40 = 0.1;            % initial vertical velocity zd (Note, x4=zd)
AoW = 15;             % initial Angle of Wings [deg]
AoV = atand(x40/x20); % initial Angle of Velocity [deg]
AoA = AoW - AoV;      % initial Angle of Attack [deg]
AoAbounds = [0 90];   % bounds for AoA [deg]


% Parameters for plots
% -------------------------------------------------------
global tt thrust drag lift grav centr aoa aow
tt = []; thrust = []; drag = []; lift = []; grav = [];
centr = []; aoa = []; aow = [];


% Call ODE solver
% -------------------------------------------------------
options = odeset('RelTol',1e-4,'AbsTol',1e-6);
[t,x] = ode23s('atmODE',[t0,tf],[x10,x20,x30,x40],options); 


% Interpolate times
% -------------------------------------------------------
% Sort times
[tt,I] = sort(tt,'ascend');
thrust = thrust(I,:);
drag = drag(I,:);
lift = lift(I,:);
grav = grav(I,:);
centr = centr(I,:);
aoa = aoa(I,:);
aow = aow(I,:);

% Average values for identical times
dt = tt(2:end)-tt(1:end-1);
I = find(dt==0);
thrust(I,:) = mean(cat(3,thrust(I,:),thrust(I+1,:)),3);
drag(I,:) = mean(cat(3,drag(I,:),drag(I+1,:)),3);
lift(I,:) = mean(cat(3,lift(I,:),lift(I+1,:)),3);
grav(I,:) = mean(cat(3,grav(I,:),grav(I+1,:)),3);
centr(I,:) = mean(cat(3,centr(I,:),centr(I+1,:)),3);
aoa(I,:) = mean(cat(3,aoa(I,:),aoa(I+1,:)),3);
aow(I,:) = mean(cat(3,aow(I,:),aow(I+1,:)),3);

J = true(length(tt),1);
J(I+1) = false;
tt = tt(J);
thrust = thrust(J,:);
drag = drag(J,:);
lift = lift(J,:);
grav = grav(J,:);
centr = centr(J,:);
aoa = aoa(J,:);
aow = aow(J,:);


% Interpolate to common times (t vector)
thrust = interp1(tt,thrust,t);
drag = interp1(tt,drag,t);
lift = interp1(tt,lift,t);
grav = interp1(tt,grav,t);
centr = interp1(tt,centr,t);
aoa = interp1(tt,aoa,t);
aow = interp1(tt,aow,t);


% Plot things
% -------------------------------------------------------
figure('Position',[0 300 600 500])
plot(x(:,3)/1000,[aow x(:,4)])
set(gca,'ylim',[0 70],'xlim',[0 80],'xtick',0:5:80,'ytick',0:5:70)
legend('AoW','Vz')
xlabel('Altitude in km')

figure('Position',[0 300 1200 800])
subplot(2,2,1)
plot(t/60,[x(:,[1 3])/1000 aow aoa]); set(gca,'ylim',[-1 100],'xlim',[t0 tf]/60)
title('positions in km'); legend('Horizontal','Vertical','AoW','AoA')
xlabel('Minutes into flight')
subplot(2,2,2)
plot(t/60,[x(:,2) x(:,4)]); %set(gca,'ylim',[-10 2500],'xlim',[t0 tf]/60)
title('velocities in m/s'); legend('Horizontal','Vertical')
xlabel('Minutes into flight')

% Separate into positive and negative contributions
thrustpos = thrust; thrustpos(thrustpos<0) = 0;
thrustneg = thrust; thrustneg(thrustneg>0) = 0;
dragpos = drag; dragpos(dragpos<0) = 0;
dragneg = drag; dragneg(dragneg>0) = 0;
liftpos = lift; liftpos(liftpos<0) = 0;
liftneg = lift; liftneg(liftneg>0) = 0;
gravpos = grav; gravpos(gravpos<0) = 0;
gravneg = grav; gravneg(gravneg>0) = 0;
centrpos = centr; centrpos(centrpos<0) = 0;
centrneg = centr; centrneg(centrneg>0) = 0;

cmap = [0 0 1; 1 0 0; 0 1 0; 0.5 0.5 0; 0 0.5 0.5];

subplot(2,2,3)
area(t,[thrustpos(:,1) dragpos(:,1) liftpos(:,1)]); hold on
h = area(t,[thrustneg(:,1) dragneg(:,1) liftneg(:,1)]); 
colormap(cmap(1:3,:)); freezeColors
h(end+1) = plot(t,thrust(:,1)+drag(:,1)+lift(:,1),'y','linewidth',2);
yl = get(gca,'ylim'); set(gca,'ylim',[max([yl(1) -20]) min([yl(2) 20])])
title('horizontal accelerations in m/s^2'); 
% legend(h,'Thrust','Drag','Lift','Sum'); 

subplot(2,2,4)
area(t,[thrustpos(:,2) dragpos(:,2) liftpos(:,2) gravpos(:,2) centrpos(:,2)]); hold on
h = area(t,[thrustneg(:,2) dragneg(:,2) liftneg(:,2) gravneg(:,2) centrneg(:,2)]);  
colormap(cmap)
h(end+1) = plot(t,thrust(:,2)+drag(:,2)+lift(:,2)+grav(:,2)+centr(:,2),'y','linewidth',2);
yl = get(gca,'ylim'); set(gca,'ylim',[max([yl(1) -20]) min([yl(2) 20])])
title('vertical accelerations in m/s^2'); 
legend(h,'Thrust','Drag','Lift','Gravity','Centrifugal','Sum','Location','NorthWest')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%