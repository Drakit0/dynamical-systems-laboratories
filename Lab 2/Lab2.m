%% Workspace initialization

clc
clear 
clear all
format short 
format compact
s = tf('s');

%% Parameters

Tsim = 0.025
L = 0.2
C = 10e-6

R = 100
% R = 500

N = 2 % Step value

I0 = 1e-2 % Initial conditions
V0 = 1
% I0 = 0
% V0 = 0

%% Responses obtained

Yf_v = tf(1,[C*L, C*R, 1])
Yl_v = tf([V0*L*C, L*I0 + V0*R*C],[C*L, C*R, 1])

Yf_i = tf([C, 0],[C*L, C*R, 1])
Yl_i = Yl_v*C*s - C*V0 

[~, impulse_len] = impulse(Yl_v, Tsim);
impulse_len = length(impulse_len) % Sampling frequency

open_system("BlockDiagram"); 
set_param("BlockDiagram", "SimulationCommand", "start"); % Execute simulation
%% Grafic resolution

tt = linspace(0, Tsim, impulse_len); % Space to plot the scope responses

figure; % v_0(t)
impulse(Yl_v, Tsim)
hold on
plot(tt, scope_response_v(:,2), DisplayName = "Scope response v(t)")
legend
hold off

figure; % I(t)
impulse(Yl_i, Tsim)
hold on
plot(tt, scope_response_i(:,2), DisplayName = "Scope response v(t)")
legend
hold off

%% Transient response verification

[n,d] = linmod("BlockDiagram");

Yf_scope_v = tf(n(1,:),d)
Yf_scope_i = tf(n(2,:),d)

transfer_functions = [Yf_v, Yf_scope_v, Yf_i, Yf_scope_i]

for i=1:1:length(transfer_functions)
    [z,p,k] = zpkdata(transfer_functions(i), "v")
    G_0 = dcgain(transfer_functions(i))
end

%% R modified

% Initialize R as 500 on the parameters section

% Execute responses obtained, graffic resolution and transient response
% verification

%% Response to a step function

% Join the step block to the block circuit as the input
% Execute responses obtained, graffic resolution and transient response
% verification

% Set the initial conditions to 0
% Execute responses obtained, graffic resolution and transient response
% verification