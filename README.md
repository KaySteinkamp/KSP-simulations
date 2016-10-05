# KSP Simulations
This is a project to design and simulate a lander using the Ion engine in Kerbal Space Program (KSP) 0.24.2.

The design part is still applicable in the newer KSP version 1.1+, but the atmospheric flight simulations assume a version less than 0.90.

# Usage
All codes are free to use for anyone.
Ion spaceship design can be optimized with `KSP_IonEngine_charts.m` and `KSP_IonSpaceship_Design.m`.
Atmospheric flight can be simulated with `KSP_AtmCalcs.m`, which contains a solver for the system of coupled non-linear ODEs describing atmospheres in KSP 0.24.2

## Ion ship design
Calculates specs for a 2-stage spaceship (traveler stage & lander stage)
powered by KSP's (0.24.2) Ion Engine.

Idea: previous studies have shown that an ion lander can actually 
achieve higher TWR (>0.5) when run on batteries rather than solar arrays,
as long as each individual burn (i.e. descend & ascend burn) lasts no 
longer than ~3 minutes. Additional benefits are that the lander can 
land at night as well as far out in the Kerbol system. It will have a
minimum amount of solar OXstats installed to recharge slowly in between
descend and ascend burns.

To collect science the lander will carry one Kerbal in an ext cmd seat as
well as lightweight science equipment. (not sure yet whether I can
include Goo or Science Jr while still achieve high enough TWR). For now I
don't plan to carry a parachute for weight reasons (only use for it would
be Duna, but I think I can make that trip without one)

The lander should be able to land (and return from) anywhere with g~<0.4, 
i.e. including Duna!
The traveling stage carries a mobile processing lab (and maybe a pod),
its main purpose is to achive high Delta-v and acceptable TWR (>0.1). It
is run on solar arrays OX4.

## Atmospheric simulations
Calculate optimal atmospheric ascend profile and/or maximum altitude and
speed obtainable with a winged spaceplane with KSP's (0.24.2) Ion Engine.
Five forces are accounted for: lift, drag, thrust, gravity, centrifugal.

1) Lift:        FL = [-vz;vy]*L*Cl*P(z)
2) Drag:        FD = -|v|*[vy;vz]*0.0049*Cd*m*P(z)
3) Thrust:      FT = T/|v|*[vy;vz]
4) Gravity:     FG = [0;-1/(R+z)^2]*G*m*M
5) Centrifugal: FC = [0;vy^2]*m/(R+z)

ODE45 is used to solve the 2 coupled (nonlinear) second-order differential
equations, by first re-writing them into a system of 4 coupled (nonlinear)
first-order differential equations.
