clear all; close all; clc;

%% Initial parameters
d = 1; % dimension
J = 1; 
h = 0; % external field
alpha = [0.1:0.1:0.9];
% ref: [Romain Bachelard, Michael Kastner]Universal Threshold for the Dynamical Behavior of Lattice Systems with Long-Range Interactions

N = round(linspace(50000,100000,4));
modo = 'static'; % random|static

time_init = 1e-9;
time_end = 1;
time_steps = 100;

repetitions = 1;

summary_mean = zeros(length(alpha),length(N),1);    
summary_std = summary_mean;

%% Create the curves and take their angular coefficients - make a fit.
for ii=1:repetitions
tic;    
	fprintf('\n%d/%d-',ii,repetitions)
    [summary_mean_temp, summary_std_temp] = generate_histogram(h,alpha,N,modo,time_init,time_end,time_steps);
    summary_mean = summary_mean + summary_mean_temp;
    summary_std = summary_std + summary_std_temp;
toc;    
end
summary_mean = summary_mean./repetitions;
summary_std = summary_std./repetitions;
%%
for ii=1:length(alpha)  
    [coeff_mean,~] = fit(log(N'),log(summary_mean(ii,:)'),'poly1');
    plot(log(N'),log(summary_mean(ii,:)')); hold on
    %[coeff_std,~]  = fit(log(N),log(summary_std(ii,:)));
    
    q_mean(ii) = coeff_mean.p1;
    %q_std(ii) = coeff_std.p1;

end





%% Plots

figure(100);
plot(alpha,q_mean,'-o')

