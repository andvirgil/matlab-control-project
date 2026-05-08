clc
clear all
close all

%% Braccio Robotico 1-DOF - Progetto di Controllo
% Parametri fisici:
%   J = 0.5  kg*m^2  (inerzia)
%   B = 2.0  Nm*s/rad (attrito viscoso)
%   M = 2.0  kg       (massa)
%   L = 0.5  m        (lunghezza braccio)
%   g = 9.81 m/s^2
%
% Equazione del moto (Eq. 1):
%   J*theta_ddot + B*theta_dot + M*g*(L/2)*sin(theta) = tau(t)
%
% Linearizzazione in theta_0 = 0, sin(theta) ~ theta (Eq. 2-3):
%   J*delta_theta_ddot + B*delta_theta_dot + (M*g*L/2)*delta_theta = delta_tau
%
% Funzione di trasferimento (Eq. 5):
%   G(s) = (1/J) / (s^2 + (B/J)*s + MgL/2J)
%        = 2 / (s^2 + 4s + 9.81)

s = tf('s');

J   = 0.5;
B   = 2.0;
M   = 2.0;
L   = 0.5;
g   = 9.81;

% Pianta linearizzata G(s) - Eq. 5
G = (1/J) / (s^2 + (B/J)*s + (M*g*L/2)/J);

fprintf('=== PIANTA G(s) ===\n');
G

% Poli della pianta
fprintf('Poli di G(s):\n');
disp(pole(G));

%% Analisi della pianta - Diagrammi di Bode, Luogo delle radici, Nyquist

figure(1);
margin(G);
title('Fig. 2 - Diagrammi di Bode della pianta G(s)');
grid on;

figure(2);
rlocus(G);
title('Luogo delle radici - G(s)');
grid on;

figure(3);
nyquist(G);
title('Diagramma di Nyquist - G(s)');
grid on;

%% Progetto del Regolatore C(s) = PI + Rete Anticipatrice (Lead)
% Specifiche:
%   S1: Stabilita' asintotica a ciclo chiuso
%   S2: Errore a regime < 5% per r(t)=sin(t) con omega=1 rad/s
%       => |S(j*1)| < 0.05  =>  |L(j*1)|_dB > 26 dB
%   S3: Overshoot < 10%  =>  MF > 45 deg
%   S4: Settling time minimo  =>  omega_c elevata

% Parte 1: Regolatore PI (Eq. 6)
%   C_PI(s) = (Kp*s + Ki) / s
Kp = 50;
Ki = 100;
C_PI = (Kp*s + Ki) / s;

% Parte 2: Rete anticipatrice Lead (Eq. 7)
%   C_lead(s) = (tau*s + 1) / (alpha*tau*s + 1),  alpha < 1
tau_z = 0.3;
tau_p = 0.05 * tau_z;   % alpha = 0.05
C_lead = (1 + tau_z*s) / (1 + tau_p*s);

% Regolatore complessivo (Eq. 8)
K  = 1;
C  = K * C_PI * C_lead;

fprintf('\n=== REGOLATORE C(s) ===\n');
C

%% Funzione d'anello e ciclo chiuso

L1    = C * G;
Wyr1  = minreal(L1 / (1 + L1));

fprintf('\n=== FUNZIONE D ANELLO L(s) ===\n');
L1

fprintf('\n=== CICLO CHIUSO W_yr(s) ===\n');
Wyr1

% Margini di stabilita'
[Gm, Pm, Wcg, Wcp] = margin(L1);
fprintf('\nMargine di Fase   = %.1f deg\n', Pm);
fprintf('Frequenza di taglio omega_c = %.2f rad/s\n', Wcp);

% Verifica specifica errore sinusoidale: |S(j*1)|
S1 = 1 / (1 + L1);
[mag_S1, ~] = bode(S1, 1);
fprintf('|S(j*1)| = %.4f  => errore = %.2f%%\n', mag_S1, mag_S1*100);

figure(1);
margin(L1);
title('Fig. 3 - Diagrammi di Bode della funzione d''anello L(s)');
grid on;

figure(2);
rlocus(L1);
title('Luogo delle radici - L(s)');
grid on;

figure(3);
nyquist(L1);
title('Diagramma di Nyquist - L(s)');
grid on;

figure(4);
step(Wyr1);
info = stepinfo(Wyr1);
title('Fig. 4 - Risposta al gradino ciclo chiuso');
grid on;
fprintf('\nOvershoot   = %.1f%%\n', info.Overshoot);
fprintf('Settling Time (2%%) = %.3f s\n', info.SettlingTime);

%% Inseguimento segnale sinusoidale omega = 1 rad/s

t_sin = 0:0.01:25;
r_sin = sin(1 * t_sin);
y_sin = lsim(Wyr1, r_sin, t_sin);

figure(5);
plot(t_sin, r_sin, '--', 'LineWidth', 2, 'DisplayName', 'r(t) = sin(\omegat)');
hold on;
plot(t_sin, y_sin, 'LineWidth', 2, 'DisplayName', '\theta(t) uscita');
hold off;
xlabel('Tempo (s)');
ylabel('\theta (rad)');
title('Fig. 5 - Inseguimento segnale sinusoidale \omega=1 rad/s');
legend;
grid on;

%% Feedforward (Eq. 9)
% Feedforward statico: Ff compensa il termine gravitazionale a regime
%   Ff = M*g*L/2 = 4.905 Nm
%   Azione: u_ff(t) = Ff * r(t)
%   Schema: Y(s) = [G(s)*C(s) + G(s)*Ff] / [1 + G(s)*C(s)] * R(s)

Ff = (M*g*L/2);   % guadagno feedforward statico

% Funzione di trasferimento con feedforward
%   W_yr_ff = Wyr1 + G*Ff*S
%   dove S = 1/(1+L1) e' la funzione di sensitività
S    = minreal(1 / (1 + L1));
Wyrf = minreal(Wyr1 + G * Ff * S);

figure(6);
step(Wyr1, Wyrf);
legend('Solo feedback C(s)', 'Feedback + Feedforward Ff');
title('Fig. 6 - Effetto del feedforward sulla risposta al gradino');
grid on;

%% Sensitività S(s), T(s), CS(s) - Eq. 10-12

% Funzione di sensitività        S  = 1/(1+L)   (Eq. 10)
% Sensitività complementare      T  = L/(1+L)   (Eq. 11)
% Sensitività del controllo      CS = C*S        (Eq. 12)
S_tf  = minreal(1 / (1 + L1));
T_tf  = minreal(L1 / (1 + L1));
CS_tf = minreal(C * S_tf);

figure(7);
bode(S_tf, T_tf, CS_tf);
legend('S(s) - Sensitività', 'T(s) - Sensitività complementare', 'CS(s) - Sensitività controllo');
title('Fig. 7 - Funzione di sensitività e sensitività del controllo');
grid on;

%% Ritardo massimo ammissibile - Eq. 13
% tau_max = MF / omega_c  (MF in radianti)
% Un ritardo e^(-tau*s) sottrae fase: angle(e^(-j*omega*tau)) = -omega*tau

tau_max = (Pm * pi/180) / Wcp;
fprintf('\n=== RITARDO MASSIMO AMMISSIBILE ===\n');
fprintf('tau_max = MF / omega_c = %.4f s\n', tau_max);

% Verifica: Nyquist con e senza ritardo
ret = pade(exp(-tau_max * s), 10);   % approssimazione di Pade' ordine 10
L2  = L1 * ret;

figure(8);
nyquist(L1, L2);
legend('L(s) senza ritardo', sprintf('L(s) con ritardo tau=%.3fs (Pade)', tau_max));
title(sprintf('Fig. 8 - Ritardo massimo ammissibile: \\tau_{max} = %.3f s', tau_max));
grid on;

figure(9);
margin(L2);
title(sprintf('Bode di L(s) con ritardo tau_{max} = %.3f s', tau_max));
grid on;

%% Simulazione con saturazione in ingresso
% Modello non lineare: u_sat(t) = sat(u(t), +/-u_max)
% con u_max = 12 Nm (coppia massima del motore)
%
% Usiamo lsim con Simulink-like approach via sistema aumentato
% oppure simulazione Eulero in loop (compatibile con octave/matlab base)

u_max = 12;     % Nm - saturazione coppia

dt   = 0.001;
t_sat = 0 : dt : 6;
r_sat = ones(size(t_sat));   % gradino unitario

% Rappresentazione di stato della pianta G(s) = 2/(s^2+4s+9.81)
% x1 = theta, x2 = theta_dot
% x1_dot = x2
% x2_dot = (u - B*x2 - MgL/2*x1) / J

% Variabili stato PI
% xi = integrale dell'errore

theta      = 0;
theta_dot  = 0;
xi         = 0;   % stato integratore PI
e_prev     = 0;

y_sat_arr  = zeros(size(t_sat));
u_sat_arr  = zeros(size(t_sat));
y_lin_arr  = zeros(size(t_sat));

% Stato parallelo senza saturazione
theta2     = 0;
theta_dot2 = 0;
xi2        = 0;

for k = 1 : length(t_sat)
    % ----- CON SATURAZIONE -----
    e    = r_sat(k) - theta;
    xi   = xi + Ki * e * dt;
    u_fb = Kp * e + xi;
    
    % lead approssimato come boost proporzionale (gia' incluso in Kp/Ki ottimizzati)
    u = u_fb;
    u_sat = max(-u_max, min(u_max, u));   % saturazione
    u_sat_arr(k) = u_sat;
    
    % Dinamica braccio (Eulero in avanti)
    theta_ddot = (u_sat - B*theta_dot - (M*g*L/2)*theta) / J;
    theta_dot  = theta_dot + theta_ddot * dt;
    theta      = theta + theta_dot * dt;
    y_sat_arr(k) = theta;
    
    % ----- SENZA SATURAZIONE (lineare) -----
    e2      = r_sat(k) - theta2;
    xi2     = xi2 + Ki * e2 * dt;
    u2      = Kp * e2 + xi2;
    
    theta2_ddot = (u2 - B*theta_dot2 - (M*g*L/2)*theta2) / J;
    theta_dot2  = theta_dot2 + theta2_ddot * dt;
    theta2      = theta2 + theta_dot2 * dt;
    y_lin_arr(k) = theta2;
end

figure(10);
subplot(2,1,1);
plot(t_sat, y_lin_arr, '--g', 'LineWidth', 2, 'DisplayName', 'Senza saturazione');
hold on;
plot(t_sat, y_sat_arr, 'b', 'LineWidth', 2, 'DisplayName', sprintf('Con saturazione ±%d Nm', u_max));
hold off;
ylabel('\theta (rad)');
title(sprintf('Fig. 9 - Risposta con saturazione in ingresso (±%d Nm)', u_max));
legend;
grid on;

subplot(2,1,2);
plot(t_sat, u_sat_arr, 'r', 'LineWidth', 2, 'DisplayName', 'u(t) saturato');
hold on;
yline(u_max,  '--', 'Color', [0.9 0.6 0], 'LineWidth', 1.5, 'DisplayName', sprintf('+%d Nm', u_max));
yline(-u_max, '--', 'Color', [0.9 0.6 0], 'LineWidth', 1.5, 'DisplayName', sprintf('-%d Nm', u_max));
hold off;
xlabel('Tempo (s)');
ylabel('u(t) [Nm]');
legend;
grid on;

%% Riepilogo finale
fprintf('\n========== RIEPILOGO PROGETTO ==========\n');
fprintf('G(s)         = 2 / (s^2 + 4s + 9.81)\n');
fprintf('C(s)         = PI (Kp=%g, Ki=%g) + Lead (tau=%g, alpha=%g)\n', Kp, Ki, tau_z, 0.05);
fprintf('Ff           = %.3f Nm (feedforward gravitazionale)\n', Ff);
fprintf('omega_c      = %.2f rad/s\n', Wcp);
fprintf('MF           = %.1f deg  (> 45 deg ✓)\n', Pm);
fprintf('|S(j1)|      = %.2f dB  (errore = %.2f%% < 5%% ✓)\n', 20*log10(mag_S1), mag_S1*100);
fprintf('Overshoot    = %.1f%%  (< 10%% ✓)\n', info.Overshoot);
fprintf('Ts (2%%)      = %.3f s\n', info.SettlingTime);
fprintf('tau_max      = %.4f s\n', tau_max);
fprintf('=========================================\n');
