function xd = atmODE(t,x)

% Define system of differential equations for y-z plane
% Note: 2nd order ODEs have been rewritten into four 1st order ODEs
% ------------------------------------------------------------------
global R M G Cd Cl_w Cl_s P0 H m T Vr AoW AoA AoV aoa aow L_w L_s AoIs
global tprev AoAbounds m0 mf mwings mctrls tburn wingDrag ctrlDrag restDrag
% Param | Value | Unit      
% R  | 6e5       | m
% M  | 5.2916e22 | kg
% G  | 6.673e-11 | grav. const.
% Cd | ?         | -
% Cl | ?         | -
% P0 | 1         | atm
% H  | 5000      | m
% m  | ?         | kg
% T  | ?         | kN
% ------------------------------------------------------------------
xd = zeros(4,1);

% Compute current mass and thrust
% -------------------------------
m = max([mf, m0 - (m0-mf)*t/tburn]);
if t>tburn
    T = 0;
end


% Terms in equations
% -------------------
AoV = atand(x(4)/x(2)); % Angle of Velocity [deg], updated in subroutine
Pair = P0*exp(-x(3)/H); % Air pressure
v = sqrt(x(2)^2+x(4)^2); % Total velocity
Grav = -G*M/(R+x(3))^2; % Gravity
Centrif = (Vr+x(2))^2/(R+x(3)); % Centrifugal acceleration

% optimize AoW/AoA within allowed range 
% (avoid negative vertical velocity and negative horizontal acceleration)
% make only small changes to AoW/AoA from previous timestep
dt = t-tprev; 
if AoA<AoAbounds(1) || AoA>AoAbounds(2)
    error('Invalid AoA from previous time step!')
end
AoAsearch = AoA + (-10:0.5:10)*(T+1)/(Cl_w+Cl_s)*dt;
AoAsearch(AoAsearch<AoAbounds(1)) = AoAbounds(1);
AoAsearch(AoAsearch>AoAbounds(2)) = AoAbounds(2);
sinAoA = sind(AoAsearch);
sinAoAi = sind(AoAsearch+AoIs);

% Compute current drag and lift coefficients
Cd = wingDrag*sinAoA/m + ctrlDrag*sinAoAi/m + restDrag*(m-mwings-mctrls)/m;
L_w = sinAoA.*(1-abs(sinAoA)).*cosd(AoAsearch);
L_s = sinAoAi;

LiftH = -x(4)*(L_w*Cl_w+L_s*Cl_s)*1e3/m*Pair; % Lift horizontal (lift-induced drag)
LiftV = x(2)*(L_w*Cl_w+L_s*Cl_s)*1e3/m*Pair;  % Lift vertical
ThrustH = T*1e3/m*cosd(AoAsearch+AoV);        % Thrust horizontal
ThrustV = T*1e3/m*sind(AoAsearch+AoV);        % Thrust vertical
DragH = -x(2)*v*0.0049*Cd*Pair;               % Drag horizontal
DragV = -x(4)*v*0.0049*Cd*Pair;               % Drag vertical
% Accelerations
xd2 = ThrustH + DragH + LiftH;
xd4 = ThrustV + DragV + LiftV + Grav + Centrif;
% Updated velocities, Delta E_kin
vvert = x(4) + xd4*dt;
%vhori = x(2) + xd2*dt;
%DEkin = (xd2.^2+xd4.^2)*dt^2 + 2*dt*(x(2)*xd2+x(4)*xd4);
% Make choice
I = find(vvert>=0 & xd2>0);
if ~isempty(I)
    % maximize sum of velocities or kinetic energy
    %[~,iopt] = max( vvert(I) + vhori(I) );
    [~,iopt] = max( 2*xd2(I) + xd4(I) );
    %[~,iopt] = max( DEkin(I) );
else
    if all(xd2<=0)
        disp('Failed to keep Ah positive')
        [~,iopt] = max( xd2(:) );
    elseif all(vvert<0)
        disp('Failed to keep Vv positive')
        [~,iopt] = max( vvert(:) );
    else
        disp('Mixed conditions')
        [~,iopt] = max( xd2(:) + xd4(:) );
    end
end
L_w = L_w(iopt);
L_s = L_s(iopt);
Cd = Cd(iopt);
LiftH = LiftH(iopt);
LiftV = LiftV(iopt);
ThrustH = ThrustH(iopt);
ThrustV = ThrustV(iopt);
DragH = DragH(iopt);
DragV = DragV(iopt);
AoA = AoAsearch(iopt);
AoW = AoA + AoV;
tprev = t;
% Differential equations
xd(1) = x(2);
xd(2) = xd2(iopt);
xd(3) = x(4);
xd(4) = xd4(iopt);


% Store parameters for plots
% ---------------------------
global thrust drag lift grav centr tt
tt = [tt;t];
thrust = [thrust;[ThrustH ThrustV]];
drag = [drag;[DragH DragV]];
lift = [lift;[LiftH LiftV]];
grav = [grav;[0 Grav]];
centr = [centr;[0 Centrif]];
aoa = [aoa;AoA];
aow = [aow;AoW];


%%% Tests %%%
disp([t tprev AoA AoW])
% disp(xd)
% if any(abs(t-(0:3600)) < 0.1)
% disp({'Time','ThrustH','ThrustV','DragH','DragV','LiftH','LiftV','Grav','Centr';...
%     t,ThrustH,ThrustV,DragH,DragV,LiftH,LiftV,Grav,Centrif})
% end
%%%%%%%%%%%%%


%----------
% dt = t-tprev;    
% AoArange = AoWrange - AoV;
% AoArange = AoArange(AoArange>AoAbounds(1) & AoArange<AoAbounds(2));
% iopt = [];
% mult = 0;
% AoAsearch = AoA;
% while dt~=0 && isempty(iopt) && all(AoAsearch>=AoArange(1)) && all(AoAsearch<=AoArange(end))
%     mult = mult + 1;
%     AoAsearch = AoA + (-1.0:0.5:1.0)*mult*dt;
%     L = sind(AoAsearch).*(1-abs(sind(AoAsearch))).*cosd(AoAsearch);
%     for i = 1:length(AoAsearch)
%         LiftH = -x(4)*L(i)*Cl*1e3/m*Pair; % Lift horizontal (lift-induced drag)
%         LiftV = x(2)*L(i)*Cl*1e3/m*Pair; % Lift vertical
%         
%         % Accelerations
%         xd2(i) = ThrustH + DragH + LiftH;
%         xd4(i) = ThrustV + DragV + LiftV + Grav + Centrif;
%         
%         % Updated velocities
%         vvert(i) = x(4) + xd4(i)*dt;
%         vhori(i) = x(2) + xd2(i)*dt;
%     end
%     iopt = find(vvert>=0 & xd2>0);
% end
% 
% if dt==0
%     L = sind(AoA).*(1-abs(sind(AoA))).*cosd(AoA);
%     LiftH = -x(4)*L*Cl*1e3/m*Pair;
%     LiftV = x(2)*L*Cl*1e3/m*Pair;
%     % Differential equations
%     xd(1) = x(2);
%     xd(2) = ThrustH + DragH + LiftH;
%     xd(3) = x(4);
%     xd(4) = ThrustV + DragV + LiftV + Grav + Centrif;
% elseif isempty(iopt)
%     disp('Failed to find AoA that fulfills conditions')
%     % guide AoA in right direction
%     if all(vvert<0)
%         % get maximum vertical lift
%         AoA = 25.7;
%     else
%         % decrease AoA to speed up horizontally
%         AoA = AoA - 10.0*dt;
%     end
%     L = sind(AoA).*(1-abs(sind(AoA))).*cosd(AoA);  
%     LiftH = -x(4)*L*Cl*1e3/m*Pair;
%     LiftV = x(2)*L*Cl*1e3/m*Pair;
%     AoW = AoA + AoV;
%     tprev = t;
%     % Differential equations
%     xd(1) = x(2);
%     xd(2) = ThrustH + DragH + LiftH;
%     xd(3) = x(4);
%     xd(4) = ThrustV + DragV + LiftV + Grav + Centrif;
% else
%     % choose AoA that maximizes Ekin
%     %[~,iopt] = max( vvert(iopt).^2 + vhori(iopt).^2 );
%     [~,iopt] = max( vvert(iopt) );
%     L = L(iopt);
%     LiftH = -x(4)*L*Cl*1e3/m*Pair;
%     LiftV = x(2)*L*Cl*1e3/m*Pair;
%     AoA = AoAsearch(iopt);
%     AoW = AoA + AoV;
%     tprev = t;
%     % Differential equations
%     xd(1) = x(2);
%     xd(2) = xd2(iopt);
%     xd(3) = x(4);
%     xd(4) = xd4(iopt);
% end
%----------