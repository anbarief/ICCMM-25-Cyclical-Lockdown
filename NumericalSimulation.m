clear all;

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
lock_down = 0;       %%% <--- SIMULATE WITH LOCKDOWN OR NOT?  1 MEANS LOCKDOWN, 0 MEANS NO LOCKDOWN
t_lag = 259;
W = 5; L = 2;        %%% <--- W DAYS OF FREE, L DAYS OF LOCKDOWN
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

function result = odefun(x, y, Damp)
  global t BETA_approx ALPHA_approx;
 [ d, ix ] = min( abs( t-x ) );
  f1 =   -(Damp(ix)*BETA_approx(ix)*y(1)*y(2));
  f2 =    Damp(ix)*(BETA_approx(ix)*y(1)*y(2)) - ALPHA_approx(ix)*y(2);
  result = [f1; f2];
end

percent_updated = 0; tic;
[t, sol] = ode45(@(t,y) odefun(t,y, Damp), t, [S(1); I(1)]);
toc;

negative = 0;
for i=1:1:nt
if sol(i,2) >= 0
S(i) = sol(i,1); I(i) = sol(i,2);
else
S(i) = S(i-1);  I(i) = 0; last_index = i; negative = 1;
break
endif
endfor
if negative == 1
for i=last_index+1:1:nt
S(i) = S(i-1);  I(i) = 0;
endfor
endif

figure;
plot(t, S, '-r', t, S_data, '--b'); legend("S Model", "S Actual");

figure;
plot(t, I, '-r', t, I_data, '--b'); legend("I Model", "I Actual");


