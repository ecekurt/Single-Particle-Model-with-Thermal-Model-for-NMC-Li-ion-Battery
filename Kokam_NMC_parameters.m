
%% Capacity

p.Ah=2.7;          %Capacity of the battery

%% Geometric Params
% Radius of particles

p.R_p = 8.5e-6;   % Radius of solid particles in positive electrode [m]
p.R_n = 12.5e-6;    % Radius of solid particles in negative electrode [m]

% Thiknesses of layers

% Fitted Kokam NMC thickness
p.L_p =8.04e-05;      % Thickness of positive electrode [m]
p.L_n =0.000182048;   % Thickness of negative electrode [m]  

%Default Kokam NMC thickness 
% p.L_p = 7e-5;      % Thickness of positive electrode [m]
% p.L_n = 7.35e-5;    % Thickness of negative electrode [m] 

% Volume fractions of electrodes

p.eps_s_p = 0.5; % Volume fraction in solid for pos. electrode
p.eps_s_n = 0.5;  % Volume fraction in solid for neg. electrode 

% Specific interfacial surface area

p.a_p = (3*p.eps_s_p)/p.R_p;  % Positive electrode [m^2/m^3]
p.a_n = (3*p.eps_s_n)/p.R_n;  % Negative electrode [m^2/m^3]

% Electrode area in C26650

p.Area_p = 0.0982;    % Positive electrode area [cm^2]
p.Area_n = 0.0982;    % Negative electrode area [cm^2]

%% Transport Params

% Diffusion coefficient in particles
% p.Ds_p0 =1e-10;  % Diffusion coeff for solid in pos. electrode, [m^2/s]
% p.Ds_n0 =1e-11;    % Diffusion coeff for solid in neg. electrode, [m^2/s] default=2e-14 

p.Ds_p0 =8e-14;  % Diffusion coeff for solid in pos. electrode, [m^2/s]
p.Ds_n0 =7e-14;    % Diffusion coeff for solid in neg. electrode, [m^2/s] 

p.Faraday = 96485.33289;      % Faraday's constant, [Coulumbs/mol]

%% Kinetic Params

p.R = 8.314472;       % Gas constant, [J/mol-K]
p.alph = 0.5;         % Charge transfer coefficients

% Reaction rates
% p.k_n0 = 3.7640e-11;  % Reaction rate in neg. electrode, [(A/m^2)*(mol^3/mol)^(1+alpha)]
% p.k_p0 = 7e-11;       % Reaction rate in pos. electrode, [(A/m^2)*(mol^3/mol)^(1+alpha)]

p.k_n0 = 1.7640e-11;  % Reaction rate in neg. electrode, [(A/m^2)*(mol^3/mol)^(1+alpha)]
p.k_p0 = 5e-11;       % Reaction rate in pos. electrode, [(A/m^2)*(mol^3/mol)^(1+alpha)]

%% Concentrations

p.c_s_p_max = 51385;   % Max concentration in cathode, [mol/m^3]  
p.c_s_n_max = 30555;    % Max concentration in anode, [mol/m^3]
p.ce = 1000;            % Fixed electrolyte concentration for SPM, [mol/m^3]

% Stochiometry

% p.theta_n_min=0.0286;              % Min soc of the negative electrode
% p.theta_n_max=0.9;            % Max soc of the negative electrode
% 
% p.theta_p_min=0.39;           % Min soc of the positive electrode
% p.theta_p_max=0.9;           % Max soc of the positive electrode

p.theta_n_min=0.0286;              % Min soc of the negative electrode
p.theta_n_max=0.4;            % Max soc of the negative electrode

p.theta_p_min=0.38;           % Min soc of the positive electrode
p.theta_p_max=0.8799;           % Max soc of the positive electrode

%% Thermodynamic Properties

p.T_amb= 298.15; 

p.Cp=750; %Heat capacity [Jkg^-1K^-1]   
p.Cc=62.7; %Heat capacity J/K (core) 
p.Cs=4.5; %Heat capacity J/K (surface)   
p.Rc=1.83;  % equiv. conduction resistance between core and surface
p.Ru=3.3;  % equiv. conduction resistance around the cell Ru=1/(h*A)
p.h=30; % 5-10 for free convection air cooling %10-70 for forced air cooling
p.A=  0.0042; %Cell area
p.M= 0.07;  % Cell mass [kg]
p.rho= 1626; 
% Heat transfer parameters
% Taken from Zhang et al (2014) [Harbin]
% http://dx.doi.org/10.1016/j.jpowsour.2014.07.110
p.C1 = 62.7;    % [J/K]
p.C2 = 4.5;     % [J/K]
p.h12 = 1.94; %1.9386; % [W/K]
p.h2a = 3.19;  % [W/K]

% Activation Energies,[J/mol]
% Taken from Prada et al (2012)
% https://doi.org/10.1149/2.064209jes
p.T_ref =298.15;       % [K]
p.E.Dsn = 35000;
p.E.Dsp =  29000;
p.E.De = 26000;
p.E.kappa_e = 34.70e3;
p.E.De_s = 37.04e3;
% Taken from Ye, Y., Shi, Y., & Tay, A. A. O. (2012).
% https://doi.org/10.1016/j.jpowsour.2012.06.055

p.E.kn= 20e3;
p.E.kp= 58e3;


%% Aging submodel parameters

p.kappa_s = 17e-5;  %1;     % [S/m] conductivity of side rxn product
p.ksei= 0.075e-14;    % [m/s^-1] neg.side kinetic rate of side rxn product adopted from Howey

p.Msei= 0.162;       % [kg/mol] molecular weight of side rxn product adopted from Delecourt
p.rhos= 100e3;        % [kg/m^3] mass density of side rxn product
p.Us= 0.4;           % [V] reference potential of side rxn
p.i0s= 1.5e-6;       % % [A/m^2] exchange current density of side rxn adopted from scott
p.L_sei=1.06e-9;        % [m] Initial thickness of the SEI layer from prada
p.Rsei_n = 1e-3;      % [Ohms*m^2] Resistivity of SEI layer, 
p.Rsei_p=0;
p.R_s_sei = 0;       % [Ohms*m^2] Resistivity of SEI layer, 
p.csolv= 4541;        % [molm^3] Solvent bulk concentration      
p.E.k_f = 6e4;        % [Jmol^-1] Activation energy of side reaction rate constant