%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Produces chart(s) for KSP's (0.24.2) Ion Engine.
%
% Motivation: I wanted to know under which conditions it is better to run
% the ion engine(s) off batteries rather than huge solar cell arrays, and
% slowly recharge batteries in between burns with very few solar cells. 
%
% Author: Kay Steinkamp, Sep 2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specs as of 0.24.2
% -------------------------------------------------------
IonEngine.thrust = 2;
IonEngine.mass = 0.25;
IonEngine.TWR = IonEngine.thrust/IonEngine.mass/9.81;
IonEngine.XpS = 0.485;
IonEngine.EpS = 8.729;
IonEngine.Isp = 4200;

Xenon.wetmass = 0.12;
Xenon.amount = 700;
Xenon.drymass = 0.05;

Battery.E = 100;
Battery.mass = 0.005;
Battery.EWR = Battery.E/Battery.mass;
Battery.Note = 'All batteries have the same Energy-to-Weight Ratio';

Solar.XL.mass = 0.35;
Solar.XL.EpS = 18;
Solar.XL.EpWS = Solar.XL.EpS/Solar.XL.mass;

Solar.OX4.mass = 0.0175;
Solar.OX4.EpS = 2;
Solar.OX4.EpWS = Solar.OX4.EpS/Solar.OX4.mass;

Solar.OXstat.mass = 0.005;
Solar.OXstat.EpS = 0.75;
Solar.OXstat.EpWS = Solar.OXstat.EpS/Solar.OXstat.mass;

Probe.mass = 0.04 + 3*0.015 + 3*0.005; 
Probe.EpS = 1.7/60;
Probe.Note = 'OKTO2; 3 LT-5 micro landing struts; 3 solar OXstat';

% Set space probe parameters
% -------------------------------------------------------
Nengine = 3;         % number of ion engines
TotBurn = 0.13*60*60;   % total burn time (sec)
IndBurn = 10:10:600; % max duration of individual burns (sec)

% Calculate parameters
% -------------------------------------------------------
IndEnergy = IndBurn.*(Probe.EpS + Nengine*IonEngine.EpS);

N.xenon = ceil(TotBurn*Nengine*IonEngine.XpS/Xenon.amount);
N.battery = ceil(IndEnergy./Battery.E);
% N.solarXL = ceil(IndEnergy./IndBurn/Solar.XL.EpS);
N.solarOX4 = ceil(IndEnergy./IndBurn/Solar.OX4.EpS);
% N.solarOXstat = ceil(IndEnergy./IndBurn/Solar.OXstat.EpS);

Mass.engine = Nengine*IonEngine.mass;
Mass.battery = N.battery*Battery.mass;
% Mass.solarXL = IndEnergy./IndBurn/Solar.XL.EpWS;
Mass.solarOX4 = N.solarOX4*Solar.OX4.mass;
% Mass.solarOXstat = IndEnergy./IndBurn/Solar.OXstat.EpWS;

Wetmass.battery = Probe.mass + Mass.engine + Mass.battery + N.xenon*Xenon.wetmass;
Wetmass.solarOX4 = Probe.mass + Mass.engine + Mass.solarOX4 + N.xenon*Xenon.wetmass;
Drymass.battery = Probe.mass + Mass.engine + Mass.battery + N.xenon*Xenon.drymass;
Drymass.solarOX4 = Probe.mass + Mass.engine + Mass.solarOX4 + N.xenon*Xenon.drymass;

TWRmin.battery = Nengine*IonEngine.thrust./Wetmass.battery/9.81;
TWRmin.solarOX4 = Nengine*IonEngine.thrust./Wetmass.solarOX4/9.81;
TWRmax.battery = Nengine*IonEngine.thrust./Drymass.battery/9.81;
TWRmax.solarOX4 = Nengine*IonEngine.thrust./Drymass.solarOX4/9.81;

Dv.battery = log(Wetmass.battery./Drymass.battery)*IonEngine.Isp*9.81;
Dv.solarOX4 = log(Wetmass.solarOX4./Drymass.solarOX4)*IonEngine.Isp*9.81;


% Chart 1 - probe parameters
% -------------------------------------------------------
figure('position',[0 0 800 500])
plot(IndBurn,repmat(N.xenon,size(IndBurn)),'y'); hold on
plot(IndBurn,N.battery,'b');
plot(IndBurn,Wetmass.battery,'b--');
plot(IndBurn,Drymass.battery,'b:');
plot(IndBurn,N.solarOX4,'r');
plot(IndBurn,Wetmass.solarOX4,'r--');
plot(IndBurn,Drymass.solarOX4,'r:');
% plot(IndBurn,N.solarOXstat,'r');
% plot(IndBurn,N.solarXL,'c');
plot(IndBurn,Dv.battery/1000,'b-.','linewidth',2)
plot(IndBurn,Dv.solarOX4/1000,'r-.','linewidth',2)
set(gca,'xtick',IndBurn([1 3:3:end]),'ylim',[0 40]);grid on
xlabel('Individual burn duration (sec)')
title(['Parameters for probe with ' num2str(Nengine) ' ion engines and ' num2str(TotBurn/60) ' minutes total burn time'])
% legend('# xenon','# batteries','# solar XL','# solar OX4','# solar OXstat')
legend('# xenon','# batteries','wet mass w/batteries','dry mass w/batteries',...
    '# solar OX4','wet mass w/solar OX4','dry mass w/solar OX4','Dv w/batteries (km/s)','Dv w/solarOX4 (km/s)')


% Chart 2 - maximum achievable TWR ("Where can we land?")
% -------------------------------------------------------
% Individual burns of 10 seconds to 10 minutes; 
% -> TWR for a probe running off batteries vs. solar (full throttle) 
% -------------------------------------------------------
figure('position',[100 100 800 500])
plot(IndBurn,TWRmin.battery,'b'); hold on
plot(IndBurn,TWRmax.battery,'b--'); hold on
% plot(IndBurn,TWRmin.solarXL,'c');
plot(IndBurn,TWRmin.solarOX4,'r');
plot(IndBurn,TWRmax.solarOX4,'r--');
% plot(IndBurn,TWRmin.solarOXstat,'r');
set(gca,'xtick',IndBurn([1 3:3:end]));grid on
xlabel('Individual burn duration (sec)')
ylabel('TWR')
title(['TWR for probe with ' num2str(Nengine) ' ion engines and ' num2str(TotBurn/60) ' minutes total burn time'])
legend('min TWR w/batteries','max TWR w/batteries','min TWR w/solar OX4','max TWR w/solar OX4')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%