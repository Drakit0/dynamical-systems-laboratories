%% Workspace initialization

clc
clear 
clear all
format short e
format compact
s = tf('s');

%% Parameter and functions definition

L = 0.2
C = 10e-6
R = 2*sqrt(L/C) % Value of R that generates two positive poles

% Pre Lab functions
% F1 = tf([1],[C*L, C*R, 1]) % V0(s)/Vi(s)
% F2 = tf([C, 0],[C*L, C*R, 1]) % I(s)/Vi(s)

[ft1_tf, ft2_tf] = RLC(R, L, C, "tf")
[ft1_zpk, ft2_zpk] = RLC(R, L, C, "zpk")
[ft1_sym, ft2_sym] = RLC(R, L, C, "sym")

ft_titles = ["ft1 tf", "ft2 tf", "ft1 zpk", "ft2 zpk", "ft1 sym", "ft2 sym"];
ft_results = [ft1_tf, ft2_tf, ft1_zpk, ft2_zpk, ft1_sym, ft2_sym];

%% Equivalence of the functions

figure;
hold on;
for i = 1:2:length(ft_results) % F1 comparation
    bode(ft_results(i));
end

grid on;
legend(ft_titles(1:2:length(ft_results)));
title("Equivalence of the 3 systems on the creation of F1")

hold off;

figure;
hold on;
for i = 2:2:length(ft_results) % F2 comparation
    bode(ft_results(i));
end

grid on;
legend(ft_titles(2:2:length(ft_results)));
title("Equivalence of the 3 systems on the creation of F2")

%% Zeros, poles, k and static gain

for i = 1:length(ft_results)
    fprintf("\nzpk and G_0 of %s \n", ft_titles{i})

    [z,p,k] = zpkdata(ft_results(i), "v")
    G_0  = dcgain(ft_results(i))
end

%% Response to a step function

N = 2;

for i = 1:length(ft_results)
    fprintf("\nresiduals and response to a step(2) function of %s \n", ft_titles{i})

    [yt, tt] = step(ft_results(i));
    yt = N*yt; % Change the step to N
    [n, d]=tfdata(ft_results(i));

    [r, p] = residue(cell2mat(n), cell2mat(d));
    figure;
    plot(tt,yt), grid, title("Reponse of " + ft_titles(i))

end

%% Modify the previous work for diff R

R = [100, 2*sqrt(L/C), 500]
r_colors = ["r", "c", "g"]

[ft_tf_r1, ~] = RLC(R(1), L, C, "tf")
[ft_tf_r2, ~] = RLC(R(2), L, C, "tf")
[ft_tf_r3, ~] = RLC(R(3), L, C, "tf")

ft_titles_r = ["ft tf R_1", "ft tf R_2", "ft tf R_3"];
ft_results_r = [ft_tf_r1, ft_tf_r2, ft_tf_r3];

figure;

for i = 1:length(ft_results_r)
    fprintf("\nzpk and G_0 of %s \n", ft_titles_r{i})

    [z,p,k] = zpkdata(ft_results_r(i), "v")
    G_0  = dcgain(ft_results_r(i))


    % fprintf("\nResponse to a step(2) function of %s \n", ft_titles{i})

    [yt, tt] = step(ft_results_r(i));
    yt = N*yt; % Change the step to N

    % [r, p] = residue(yt, tt);
    hold on
    plot(tt,yt, color = r_colors(i))
end

legend(ft_titles_r);
grid on;
title("Reponse with differen R to a step of value 2");

%% Free response

I_0 = 5e-3
V_0 = 2

for i = 1:length(R)
    fprintf("\nFree response with %d ohms", R(i))
    Yl = tf(L*[C*V_0, I_0] + C*R(i)*V_0,[C*L, C*R(i), 1])
end

%% Resolution of 2.7

% EDO i(t) + C(dv(t)/dt) = vi(t)/R - i(t) - L/R * (di(t)/dt)

R = 0.524 % Greatest R value
L = 0.2
C = 5e-6

V_0 = 2
I_0 = 5/(2*R)

N = 5

Yf = tf(1/R, [1, 1/(5*R) + 5*R, 2])
[yfu,tt] = step(Yf);
yfu = yfu*N;

Yl = tf([1/R, (1/R)*(1/(5*R) + 5*R)], [1, 1/(5*R) + 5*R, 2])
yl_t = impulse(Yl);


plot(tt, yfu + yl_t);
grid on;
title("y_s(t) of the problem 2.7 for the greatest R")
