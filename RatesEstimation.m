%%Arief Anbiya (15-2-2025)

clear all;

write_to_csv = 1;

%IMPORTING DATA
S_approx = csvread("S_ready_USA.csv", 1, 0);
I_approx = csvread("I_ready_USA.csv", 1, 0);
Mu = csvread("Mu_ready_USA.csv", 1, 0);
t_samples = csvread("index_ready_USA.csv", 1, 0);
%%%%%%%%%%%%%%%

dt=1/100; nt = length(t_samples);
time_index = transpose([0:0.25:t_samples(end)]); NT = length(time_index);
locator = transpose([1:1:nt]);

Sigma = Mu.*S_approx;
Lambda = Mu.*I_approx;
N_p = S_approx(1);

%%% Normalize the data
factor = 1/N_p;
S_initial = factor*S_approx;
I_initial = factor*I_approx;
Sigma = factor*Sigma;
Lambda = factor*Lambda;
%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
err_tol = 10^(-6);
err_tol2 = 10^(-9);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

percent_updated = 0;
time_start = 1;

for k = [1:1:nt]

if 1 + k*4 > NT
  disp("Finished");
  nb = length(BETA);
  remaining = nt - nb;
  BETA(nb+1:nb+remaining) = BETA(nb);
  GAMMA(nb+1:nb+remaining) = GAMMA(nb);
  break
else
  start_index = time_start;
  end_index = sum(locator.*(time_index(1 + k*4)==t_samples ));
  S_guess = S_initial(start_index:end_index);
  I_guess = I_initial(start_index:end_index);
  mu_sub = Mu(start_index:end_index);
  sigma_sub = Sigma(start_index:end_index);
  lambda_sub = Lambda(start_index:end_index);

  beta_guess = 0.001;
  gamma_guess = 0.001;

  n_mat = length(S_guess);

  if k == 1
    M = zeros(n_mat, n_mat);
    M(1,1) = -1;
    M(n_mat, n_mat) = 1;
    for i=[2:(n_mat-1)]
      M(i,i) = -(2 + dt*mu_sub(i));
    endfor

    for i=[1:1:n_mat-1]
    M(i, 1+i) = 1;
    M(1+i, i) = 1;
    endfor

    M(n_mat, n_mat-1) = -1;
    INV = inv(M);
  endif

  error_beta = 1; error_gamma = 1;

while ((error_beta > err_tol2) || (error_gamma > err_tol2))

  error_S = 1; error_I = 1;

while ((error_S > err_tol) || (error_I > err_tol))
S_RHS(1,1) = -beta_guess*(  (0.5*(S_guess(2)+S_guess(1)))*(0.5*(I_guess(2)+I_guess(1))));
S_RHS(2:n_mat-1,1) = -beta_guess*(( (0.5*(S_guess(3:n_mat)+S_guess(2:n_mat-1))).*(0.5*(I_guess(3:n_mat)+I_guess(2:n_mat-1)))) ...
- (0.5*(S_guess(1:n_mat-2)+S_guess(2:n_mat-1))).*(0.5*(I_guess(1:n_mat-2)+I_guess(2:n_mat-1)))) ...
- (mu_sub(2:n_mat-1).*sigma_sub(2:n_mat-1))  ;
S_RHS(n_mat,1) = -beta_guess*(  (0.5*(S_guess(n_mat)+S_guess(n_mat-1)))*(0.5*(I_guess(n_mat)+I_guess(n_mat-1))));
new_S = INV*S_RHS*dt;

I_RHS(1,1) = beta_guess*(  (0.5*(S_guess(2)+S_guess(1)))*(0.5*(I_guess(2)+I_guess(1)))) - gamma_guess*(0.5*(I_guess(2)+I_guess(1)));
I_RHS(2:n_mat-1,1) = beta_guess*(( (0.5*(S_guess(3:n_mat)+S_guess(2:n_mat-1))).*(0.5*(I_guess(3:n_mat)+I_guess(2:n_mat-1)))) ...
- (0.5*(S_guess(1:n_mat-2)+S_guess(2:n_mat-1))).*(0.5*(I_guess(1:n_mat-2)+I_guess(2:n_mat-1)))) ...
- gamma_guess*(0.5*(I_guess(3:n_mat)+I_guess(2:n_mat-1)) - 0.5*(I_guess(1:n_mat-2)+I_guess(2:n_mat-1))) ...
- mu_sub(2:n_mat-1).*lambda_sub(2:n_mat-1);
I_RHS(n_mat,1) = beta_guess*( (0.5*(S_guess(n_mat)+S_guess(n_mat-1)))*(0.5*(I_guess(n_mat)+I_guess(n_mat-1)))) - gamma_guess*(0.5*(I_guess(n_mat)+I_guess(n_mat-1)));
new_I = INV*I_RHS*dt;

error_S = sqrt(sum((S_guess - new_S).^2))/sqrt(sum(new_S.^2));
error_I = sqrt(sum((I_guess - new_I).^2))/sqrt(sum(new_I.^2));

S_guess = new_S;
I_guess = new_I;

end

if sum(I_guess == 0) > 1
  beta_new = 0;
  gamma_new = 0;
else
alpha_10 = -2*sum( (0.25)*(I_guess(1:n_mat-1)+I_guess(2:n_mat)).*(S_guess(1:n_mat-1)+S_guess(2:n_mat)).* ...
 (I_guess(2:n_mat)-I_guess(1:n_mat-1) - (S_guess(2:n_mat)-S_guess(1:n_mat-1)) ));
alpha_01 = 2*sum( 0.5*(I_guess(1:n_mat-1)+I_guess(2:n_mat)).*(I_guess(2:n_mat)-I_guess(1:n_mat-1)) );

alpha_20 = 2*(dt)*sum(  (0.25^2)*((I_guess(1:n_mat-1)+I_guess(2:n_mat)).*(S_guess(1:n_mat-1)+S_guess(2:n_mat))).^2 );
alpha_11 = -2*(dt)*sum( (0.25*(I_guess(1:n_mat-1)+I_guess(2:n_mat)).^2).*(0.5*(S_guess(1:n_mat-1)+S_guess(2:n_mat)))  );
alpha_02 = (dt)*sum( (0.25*(I_guess(1:n_mat-1)+I_guess(2:n_mat)).^2) );

beta_new =-((2*alpha_02*alpha_10) - (alpha_01*alpha_11))/(-(alpha_11^2) + 4*alpha_02*alpha_20);
gamma_new =-((alpha_10*alpha_11) - 2*(alpha_01*alpha_20))/((alpha_11^2) - 4*alpha_02*alpha_20);
end

error_beta = abs(beta_new - beta_guess);
error_gamma = abs(gamma_new - gamma_guess);

beta_guess = beta_new;
gamma_guess = gamma_new;
end

BETA_daily(k) = beta_new;
GAMMA_daily(k) = gamma_new;
BETA(start_index:end_index) = beta_new;
GAMMA(start_index:end_index) = gamma_new;

percent = floor(100*end_index/nt);
if ((mod(percent,5)==0) && (percent ~= percent_updated))
      disp([num2str(percent), '%, beta = ', num2str(beta_new)]);
      percent_updated = percent;
endif

time_start = end_index;
end
endfor

n_day = length(BETA_daily);
t_day = [0:1:n_day];

BETA_daily = [BETA_daily BETA_daily(end)];
GAMMA_daily = [GAMMA_daily GAMMA_daily(end)];
n_daily = length(BETA_daily);

BETA_interp = interp1(t_day, BETA_daily, t_samples);
GAMMA_interp = interp1(t_day, GAMMA_daily, t_samples);

figure;
plot(t_samples, BETA_interp); title("Beta*Np");
figure;
plot(t_samples, GAMMA_interp); title("Gamma (or Alpha)");

if write_to_csv == 1
csvwrite("Beta.csv", BETA_interp/N_pop);
csvwrite("Gamma.csv", GAMMA_interp);
endif
