clear all
close all
clc

%% Simulation Inputs and Initial Conditions


% Define Constant Values
g_k = 36; %[mS/cm^2]
g_na = 120; %[mS/cm^2]
g_l = 0.3; %[mS/cm^2]

C_m = 1; %[uF/cm^2]

E_k = -12; %[mV]
E_na = 115; %[mV]
E_l = 10.6; %[mV]

Vrest = -70; %[mV]
VV = 0; % Simulation Voltage which will be corrected after calcuations

%Define Time Step
ts = .01; %[s]
tt = 0:ts:100; 


%Define Input Current
II = zeros(1,length(tt));
II(1:5/ts) = 5; %Vector defines input I for each time step
%% Simulation Iteration

% Initialize Variables with Steady State Values
alpha_m = .1*((25-VV)/(exp((25-VV)/10)-1));
beta_m = 4*exp(-VV/18);
alpha_n = .01*((10-VV)/(exp((10-VV)/10)-1));
beta_n = .125*exp(-VV/80);
alpha_h = .07*exp(-VV/20);
beta_h = 1/(exp((30-VV)/10)+1);

% Calculate Voltage Gated Channel Probabilities
nn = alpha_n/(alpha_n+beta_n); 
mm = alpha_m/(alpha_m+beta_m);
hh = alpha_h/(alpha_h+beta_h);

nvec = nn;
mvec = mm;
hvec = hh;
Vol = [VV+Vrest];

% Calcuate Initial Conductances
g_na_vec = mm^3*g_na*hh;
g_k_vec = nn^4*g_k;

% Iterate over time
for ii = 1:length(tt)-1
    
    
    alpha_m = .1*((25-VV)/(exp((25-VV)/10)-1));
    beta_m = 4*exp(-VV/18);
    alpha_n = .01*((10-VV)/(exp((10-VV)/10)-1));
    beta_n = .125*exp(-VV/80);
    alpha_h = .07*exp(-VV/20);
    beta_h = 1/(exp((30-VV)/10)+1);
    
    % Current time derivative for each probability
    dmm = alpha_m*(1-mm)-beta_m*mm;
    dnn = alpha_n*(1-nn)-beta_n*nn;
    dhh = alpha_h*(1-hh)-beta_h*hh;
    
    %Update conductances over time
    g_na_vec = [g_na_vec, mm^3*g_na*hh];
    g_k_vec = [g_k_vec, nn^4*g_k]; 
    
    
    nvec = [nvec, nn+ts*dnn];
    nn  = nn+ts*dnn;
    mvec  = [mvec, mm+ts*dmm];
    mm = mm+ts*dmm;
    hvec  = [hvec, hh+ts*dhh];
    hh = hh+ts*dhh;
    
    % Calculate Ion Currents
    I_na = mm^3*g_na*hh*(VV-E_na);
    I_k = nn^4*g_k*(VV-E_k);
    I_l = g_l*(VV-E_l);
    
    I_ion = II(ii) - I_k - I_na - I_l;
   
    % Current time derivative of membrane motental
    dVV = I_ion/C_m;
    
    % Update membrane potential
    Vol = [Vol, (VV+ts*dVV) + Vrest];
    VV = (VV+ts*dVV);
    
end

%% Plot voltage results

figure(1)
plot(tt,Vol)
xlabel('Time [ms]')
ylabel('Membrane Potential [mV]')

figure(2)
plot(tt, g_na_vec, tt, g_k_vec)
legend('g_Na','g_K')
xlabel('Time [ms]')
ylabel('Conductance [mS/cm^2]')
    



