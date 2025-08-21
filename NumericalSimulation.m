%Author: Arief Anbiya

clear all; close all;

global t BETA_approx ALPHA_approx;

dt = 1/100;
t = csvread("index_ready_USA.csv", 1, 0);
nt = length(t);
time_index = transpose([0:0.25:t(end)]);
I_data = csvread("I_ready_USA.csv", 1, 0);
S_data = csvread("S_ready_USA.csv", 1, 0);
BETA = csvread("Beta.csv");
ALPHA = csvread("Alpha.csv");

%USA 1 March 2020
S = zeros(nt,1); I = zeros(nt,1);
S(1) = S_data(1); I(1) = I_data(1);

BETA_approx = BETA; ALPHA_approx = ALPHA;

%%%%%%%%%%%%%% LOCK DOWN PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% EDIT THESE 4 LINES FOR DIFFERENT SCENARIOS:
lock_down = 1;       %%% <--- SIMULATE WITH LOCKDOWN OR NOT?  1 MEANS CYCLICAL LOCKDOWN, 0 MEANS NO LOCKDOWN WHICH IS FITTING OF THE ACTUAL DATA
t_lag = 259;
W = 4; L = 3;        %%% <--- W DAYS OF FREE, L DAYS OF LOCKDOWN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% Create the Damping vector
checkpoint = t_lag;
for it = 1:1:nt
if t(it) <= t_lag
Damp(it) = 1;
else
  if t(it) >= checkpoint + W+L
    checkpoint = checkpoint + W+L;
  endif
  if (checkpoint <= t(it)) &&  (t(it) <= checkpoint + W)
    Damp(it) = 1;
  else
    Damp(it) = 0.175;
  endif
endif
endfor

if lock_down == 0
Damp = ones(nt,1);
end
figure;
plot(t, Damp); title("Damping function");
ylim([-0.1, 1.1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Runge-Kutta 4 implementation
function result = f(y, x, Damp)
  global t BETA_approx ALPHA_approx;
  [ d, ix ] = min( abs( t-x ) );
  f1 =   -(Damp(ix)*BETA_approx(ix)*y(1)*y(2));
  f2 =    Damp(ix)*(BETA_approx(ix)*y(1)*y(2)) - ALPHA_approx(ix)*y(2);
  result = [f1; f2];
end
percent_updated = 0;
tic;
for it=[2:1:nt]
yn = [S(it-1); I(it-1)];
k1 = dt*f( yn, t(it-1), Damp );
k2 = dt*f( yn + 0.5*k1, t(it-1) + 0.5*dt, Damp );
k3 = dt*f( yn + 0.5*k2, t(it-1) + 0.5*dt, Damp );
k4 = dt*f( yn + k3, t(it-1) + dt, Damp  );
yn_new = yn + (k1/6) + (k2/3) + (k3/3) + (k4/6);
S(it) = yn_new(1); I(it) = yn_new(2);
if I(it) <= 0
I(it) = 0;
end
percent = floor(100*it/nt);
if ((mod(percent,5)==0) && (percent ~= percent_updated))
      disp(['Computing... ', num2str(percent), '%']);
      percent_updated = percent;
endif
endfor
toc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

figure;
plot(t, S, '-r', t, S_data, '--b'); legend("S Model", "S Actual");

figure;
plot(t, I, '-r', t, I_data, '--b'); legend("I Model", "I Actual");


