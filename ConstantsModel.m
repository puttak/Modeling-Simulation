clc
clear all
close all

% Defining the require constans of catalyts and reactor dimensions

L = 2.5;                     % [m]              Length of reactor
dt = 0.0256;                 % [m]              Diameter of reactor
dp = 0.0082;                 % [m]              Diameter of catalyts particle
epsilon = 0.48;              % [m^3/m^3]        Porosity
Density_bed = 50;            % [kg/m^3]         Density of reactor bed

% Defining the require constans of operation conditions

P = 1;                        % [atm]           Pressure of reactor bed
Tb = 450:600;                 % [oC]            Temperature of reactor coolant
T0 = 300;                     % [oC]            Temperature of inlet reactor
Rep = 1400;                   % [-]             Reynolds number
Flowin = 4;                   % [Nm^3/h]        Inlet volume flowrate 
C_ethan = 1:2;                % [%mol]          Mole frac of inlet ethane
c_air = 98:99;                % [%mol]          Mole frac of inlet Air

% Defining the require constans of transfer parameters

Deffr = 32;                   % [m^2/h]         Effective mass transfer 
                              %                 coefficient in radius direction
Deffz = 53;                   % [m^2/h]         Effective mass transfer coefficient 
                              %                 in horizontal axis direction
kg = 576;                     % [m^3/(m^2*h)]   Surface mass transfer coefficient
hg = 928.8;                   % [kJ/(m^2*h*K)]  Surface heat transfer coefficient
keffr = 9.72;                 % [kJ/(m*h*K)]    Effective thermal conductivity 
                              %                 in the radius direction
hw = 1051.2;                  % [kJ/(m^2*h*K)]  Wall heat transfer coefficient

% Defining the require constants properties of components

C2H6 = struct('Mw',30.07,   'Tc',305.406,   'Pc',4880109,  'cp',[1.131,0.019225,-0.000005561,0,1500]);
C2H4 = struct('Mw',28.054,  'Tc',282.3438,  'Pc',5045427,  'cp',[1.424,0.014394,-0.000004392,0,1500]);
O2   = struct('Mw',31.998,  'Tc',154.645,   'Pc',5043213,  'cp',[3.639,0.000506,0,-22700,2000]);
CO2  = struct('Mw',44.009,  'Tc',304.1548,  'Pc',7380862,  'cp',[5.457,0.001045,0,-115700,2000]);
CO   = struct('Mw',28.01,   'Tc',134.18,    'Pc',3710046,  'cp',[3.376,0.000557,0,-3100,2500]);
H2O  = struct('Mw',18.015,  'Tc',647.1081,  'Pc',22072227, 'cp',[3.47,0.00145,0,12100,2000]);

components = [C2H6 C2H4 O2 CO2 CO H2O];




