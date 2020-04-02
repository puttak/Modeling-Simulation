clc;clear all;close all;
%% Defining the require constans of catalyts and reactor dimensions

L = 2.5;                         % [m]                              Length of reactor
dt = 0.0256;                     % [m]                              Diameter of reactor
dp = 0.0082;                     % [m]                              Diameter of catalyts particle
epsilon = 0.48;                  % [m^3/m^3]                        Porosity
Density_bed = 50*(1000);         % [kg/m^3][g/m^3]                  Density of reactor bed
Phis = 1;                        % [unitless]                       Sphericity factor for spheres catalyts particles 
as = 6*(1-epsilon)/(Phis*dp);    % [m^2/m^3]                        External surface to particle volume ratio

%% Defining the require constans of operation conditions

P = 1*(101325);                  % [atm][Pa]                        Pressure of reactor bed
Tb = 450+(273.15):600+(273.15);  % [oC][K]                          Temperature of reactor coolant
T0 = 300+(273.15);               % [oC][K]                          Temperature of inlet reactor
Rep = 1400;                      % [unitless]                       Reynolds number
Flowin = 4*(1/3600);             % [Nm^3/h][Nm^3/s]                 Inlet volume flowrate 
C_ethan = 1:2;                   % [%mol]                           Mole frac of inlet ethane
c_air = 98:99;                   % [%mol]                           Mole frac of inlet Air

%% Defining the require constans of transfer parameters

Deffr = 32*(1/3600);             % [m^2/h][m^2/s]                   Effective mass transfer 
                                 %                                  coefficient in radius direction.
Deffz = 53*(1/3600);             % [m^2/h][m^2/s]                   Effective mass transfer coefficient 
                                 %                                  in horizontal axis direction.
kg = 576*(1/3600);               % [m^3/(m^2*h)][m^3/(m^2*s)]       Surface mass transfer coefficient
hg = 928.8*(1000)*(1/3600);      % [kJ/(m^2*h*K)][J/(m^2*s*K)]      Surface heat transfer coefficient
keffr = 9.72*(1000)*(1/3600);    % [kJ/(m*h*K)][J/(m*s*K)]          Effective thermal conductivity 
                                 %                                  in the radius direction.
hw = 1051.2*(1000)*(1/3600);     % [kJ/(m^2*h*K)][J/(m^2*s*K)]      Wall heat transfer coefficient

%% Defining the require constants of components properties 
%
% Units
%
% Mw: [g/mol]           Tc: [K]       Pc: [Pa]      
% cp_R: [unitless depend on R]        deltaS0: [J/(mol*K)]      deltaH0:[kJ/mol][J/mol]
% R=8.314 (J/mol*K)

C2H6 = struct('Mw',30.07,   'Tc',305.406,   'Pc',4880109,...
    'cp_R',[1.131,0.019225,-0.000005561,0,1500] ,'deltaS0',5.27e01 ,'deltaH0',(1000)*4.8e01);
C2H4 = struct('Mw',28.054,  'Tc',282.3438,  'Pc',5045427,...
    'cp_R',[1.424,0.014394,-0.000004392,0,1500] ,'deltaS0',4.34e01 ,'deltaH0',(1000)*1.48e02);
O2   = struct('Mw',31.998,  'Tc',154.645,   'Pc',5043213,...
    'cp_R',[3.639,0.000506,0,-22700,2000]       ,'deltaS0',5.59e01 ,'deltaH0',(1000)*6.02e01);
CO2  = struct('Mw',44.009,  'Tc',304.1548,  'Pc',7380862,...
    'cp_R',[5.457,0.001045,0,-115700,2000]      ,'deltaS0',5.66e01 ,'deltaH0',(1000)*8.38e01);
CO   = struct('Mw',28.01,   'Tc',134.18,    'Pc',3710046,...
    'cp_R',[3.376,0.000557,0,-3100,2500]        ,'deltaS0',8.66e01 ,'deltaH0',(1000)*4.09e01);
H2O  = struct('Mw',18.015,  'Tc',647.1081,  'Pc',22072227,...
    'cp_R',[3.47,0.00145,0,12100,2000]          ,'deltaS0',5.27e01 ,'deltaH0',(1000)*8.63e01);
components = [C2H6 C2H4 O2 CO2 CO H2O];

%% Defining the require constants for kinetic of reactions
%
% Units
%
% Aprime(A'): [mmol/(g*h)][mol/(g*s)]   EnergyA: [kJ/mol][J/mol]      m: [unitless]
% component coefficients  vcoffrxn: [unitless]
% component order list: [C2H6 C2H4 O2 CO2 CO H2O]

RnxKinetic = struct('Aprime',[4.95 1.35 1.76 2.61 2.16]*(1/1000)*(1/3600),...
    'EnergyA', [7.55e01 5.24e01 1.43e02 1.10e02 8.80e01]*(1000),...  
    'm', [1 5.45e-02 1.07 1.71e-01 5.38e-01],...
    'vcoffrxn', [-1 1 -0.5 0 0 1; -1 0 -3.5 2 0 3; -1 0 -2.5 0 2 3; 0 -1 -3 2 0 2; 0 -1 -2 0 2 2]); 

%% Defining the require constants for 

%% Calculation 
Nz = 5; % No. of interior point in z direction
Nr = 3; % No. of interior point in r direction
zmin = 0; zmax = 1;
rmin = 0; rmax = 1;
z_nodes = [0,0.001,sort(Roots_of_Jacobi_Polynomial(0,0,Nz))',0.999,1] ;  % Roots of Jacobi polynomial with (a,b==0) in z direction
z_nodes = (zmax-zmin)*z_nodes+zmin;
r_nodes = [0,0.001,sort(Roots_of_Jacobi_Polynomial(0,0,Nr))',0.999,1] ;  % Roots of Jacobi polynomial with (a,b==0) in r direction
r_nodes = (rmax-rmin)*r_nodes+rmin;
syms z
Lz = sym(ones(numel(z_nodes),1));
for i=1:numel(z_nodes)
    for j=1:numel(z_nodes)
        if j~=i
            Lz(i,1) = (z-z_nodes(j))/(z_nodes(i)-z_nodes(j))*Lz(i,1);
            % Lz is Lagrange interpolation polynomial in z direction
        end
    end 
end
syms r
Lr = sym(ones(numel(r_nodes),1));
for i=1:numel(r_nodes)
    for j=1:numel(r_nodes)
        if j~=i
            Lr(i,1) = (r-r_nodes(j))/(r_nodes(i)-r_nodes(j))*Lr(i,1);
            % Lr is Lagrange interpolation polynomial in r direction
        end
    end 
end
Lz_prime = diff(Lz);      % First drivative of Lagrange polynomial in z direction
Lz_Zegond = diff(Lz,2);   % Second drivative of Lagrange polynomial in z direction
Lr_prime = diff(Lr);      % First drivative of Lagrange polynomial in r direction
Lr_Zegond = diff(Lr,2);   % Second drivative of Lagrange polynomial in r direction
Az = zeros(numel(z_nodes));
Bz = zeros(numel(z_nodes));
for i = 1:numel(z_nodes)
    for j = 1:numel(z_nodes)
        Az(i,j) = double(subs(Lz_prime(j),z_nodes(i)));
        Bz(i,j) = double(subs(Lz_Zegond(j),z_nodes(i))); 
    end    
end
Ar = zeros(numel(r_nodes));
Br = zeros(numel(r_nodes));
for i = 1:numel(r_nodes)
    for j = 1:numel(r_nodes)
        Ar(i,j) = double(subs(Lr_prime(j),r_nodes(i)));
        Br(i,j) = double(subs(Lr_Zegond(j),r_nodes(i)));
    end    
end

















