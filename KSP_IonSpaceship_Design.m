%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculate specs for a 2-stage spaceship (traveler stage & lander stage)
% powered by KSP's (0.24.2) Ion Engine.
%
% Idea: previous studies have shown that an ion lander can actually 
% achieve higher TWR (>0.5) when run on batteries rather than solar arrays,
% as long as each individual burn (i.e. descend & ascend burn) lasts no 
% longer than ~3 minutes. Additional benefits are that the lander can 
% land at night as well as far out in the Kerbol system. It will have a
% minimum amount of solar OXstats installed to recharge slowly in between
% descend and ascend burns.
%
% To collect science the lander will carry one Kerbal in an ext cmd seat as
% well as lightweight science equipment. (not sure yet whether I can
% include Goo or Science Jr while still achieve high enough TWR). For now I
% don't plan to carry a parachute for weight reasons (only use for it would
% be Duna, but I think I can make that trip without one)
%
% The lander should be able to land (and return from) anywhere with g~<0.4, 
% i.e. including Duna!
% The traveling stage carries a mobile processing lab (and maybe a pod),
% its main purpose is to achive high Delta-v and acceptable TWR (>0.1). It
% is run on solar arrays OX4.
%
% Author: Kay Steinkamp, Sep 2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Part Specs as of 0.24.2 
% (ignore massless parts - wannabe realistic)
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


%% Set fixed parameters
% -------------------------------------------------------
% ScanSAT adapter
sat.parts = {
    '4 solar OXstat';
    '2 Jr. Docking Ports';
    '1 Radar low-res';
    '1 Multispec';
    '1 Karbonite'};
sat.mass = sum([...
    4*0.005;
    2*0.02;
    0.03;
    0.03;
    0.05]);
    
% Science adapter
sci.parts = {
    '2 Jr. Docking Ports';
    '1 Science Jr.';
    '2 Goo';
    '1 Sensor Nose Cone';
    '2 Illuminators'};
sci.mass = sum([...
    2*0.02;
    0.2;
    2*0.15;
    0.08;
    2*0.02]);

%% Landing Stage
LS.Nengine = 3;     % number of ion engines
LS.totBurn = 8*60;  % total burn time (sec) for descent+ascend
LS.indBurn = 0.67*LS.totBurn;  % max duration of individual burns (sec)
                               % (assume ascend takes 2/3 of total burn
                               % time for Duna (atmosphere!))
LS.fixedParts = {
    '1 OKTO 2'; 
    '3 LT-5 micro landing struts'; 
    '3 solar OXstat (M=0?)';
    '1 ext cmd seat';
    '1 Kerbal';
    '1 PresMat barometer (M=0?)';
    '1 GRAVMAX detector (M=0?)';
    '1 Seismic accelerometer (M=0?)';
    '1 2HOT thermometer (M=0?)';
    '1 Communotron 16 (M=0?)';
    '1 Jr. Docking Port'};
LS.fixedMass = sum([...
    0.04;
    3*0.015;
    3*0.005;
    0.05;
    0.09375;
    0.005;
    0.005;
    0.005;
    0.005;
    0.005;
    0.02]);


%% Traveling Stage
TS.Nengine = 8;
TS.totBurn = 6*60*60; 
TS.fixedParts = {
    '1 OKTO 2';
    '1 aSAS';
    '1 Communotron';
    '4 Illuminators';
    'Batteries (400)';
    '1 Docking Port';
    '2 Jr. Docking Ports';
    'KAS container (filled ?0.2t)';
    'RCS fuel (200)';
    '4 RCS thrusters';
    '1 MPL (Mobile Processing Lab)';
    'MPL: 4 parachutes';
    'MPL: 4 LT-5 micro landing struts'
    'MPL: 1 OKTO 2';
    'MPL: 4 RCS thrusters';
    'MPL: batteries (200)';
    'MPL: RCS fuel (50)';
    'MPL: 1 Long Ladder';
    'MPL: 1 Docking Port'};
TS.fixedMass = sum([...
    0.04;
    0.1;
    0.005;
    4*0.02;
    0.02;
    0.05;
    2*0.02;
    0.2;
    1.1;
    4*0.05;
    3.5;
    4*0.15;
    4*0.015;
    0.04;
    4*0.05;
    0.01;
    0.25;
    0.005;
    0.05]);
    

%% Calculate parameters
% -------------------------------------------------------
% Landing Stage
LS.indEnergy = LS.indBurn*LS.Nengine*IonEngine.EpS;
LS.Nxenon = ceil(LS.totBurn*LS.Nengine*IonEngine.XpS/Xenon.amount);
LS.Nbattery = ceil(LS.indEnergy/Battery.E);

LS.wetMass = LS.Nengine*IonEngine.mass + LS.Nbattery*Battery.mass + ...
    LS.Nxenon*Xenon.wetmass + LS.fixedMass;
LS.dryMass = LS.Nengine*IonEngine.mass + LS.Nbattery*Battery.mass + ...
    LS.Nxenon*Xenon.drymass + LS.fixedMass;

LS.TWR = LS.Nengine*IonEngine.thrust/9.81 ./ [LS.wetMass LS.dryMass];

LS.Dv = log(LS.wetMass/LS.dryMass)*IonEngine.Isp*9.81;


% Traveling Stage
TS.NsolarOX4 = ceil((TS.Nengine+LS.Nengine)*IonEngine.EpS/Solar.OX4.EpS);
TS.Nxenon = ceil(TS.totBurn*(TS.Nengine+LS.Nengine)*IonEngine.XpS/Xenon.amount);

TS.wetMass = TS.Nengine*IonEngine.mass + TS.NsolarOX4*Solar.OX4.mass + ...
    TS.Nxenon*Xenon.wetmass + TS.fixedMass + LS.wetMass + sat.mass + sci.mass;
TS.dryMass = TS.Nengine*IonEngine.mass + TS.NsolarOX4*Solar.OX4.mass + ...
    TS.Nxenon*Xenon.drymass + TS.fixedMass + LS.dryMass + sat.mass + sci.mass;

TS.TWR = (TS.Nengine+LS.Nengine)*IonEngine.thrust/9.81 ./ [TS.wetMass TS.dryMass];

TS.Dv = log(TS.wetMass/TS.dryMass)*IonEngine.Isp*9.81;


disp(sat)
disp(sci)
disp(LS)
disp(TS)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%