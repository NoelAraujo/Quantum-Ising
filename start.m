clear all; close all; clc;

%% Initial parameters
d = 1; % dimension
J = 1; 
h = 0; % external field
alpha = [0.4 0.8];
% ref: [Romain Bachelard, Michael Kastner]Universal Threshold for the Dynamical Behavior of Lattice Systems with Long-Range Interactions

N = [100 200 300];
modo = 'random'; % random|static

time_init = 1e-9;
time_end = 1;
time_steps = 250;

repetitions = 5;

summary_mean = zeros(length(alpha),length(N),1);    
summary_std = summary_mean;
tic;
for ii=1:repetitions
    
	fprintf('\n%d/%d-',ii,repetitions)
    [summary_mean_temp, summary_std_temp] = generate_histogram(h,alpha,N,modo,time_init,time_end,time_steps);
    summary_mean = summary_mean + summary_mean_temp;
    summary_std = summary_std + summary_std_temp;
    
end
toc;
summary_mean = summary_mean./repetitions;
summary_std = summary_std./repetitions;



%% Plots
pack_colors = ['r','b','k','g'];


close all
figure(1);
show_this(pack_colors, 'o-', alpha, log(N),log(summary_mean),'mean','log(N)','<t_0>')

figure(2);
show_this(pack_colors, 'o-',alpha, log(N),log(summary_std),'sigma','log(N)','std')

figure(3);
for ii=1:length(alpha)
    h = plot(log(N),log(summary_std(ii,:,1)./summary_mean(ii,:,1)),strcat(pack_colors(ii),'o-'));
    set(h, 'MarkerFaceColor', get(h, 'Color'));
    hold on
    Legenda{ii} = strcat('\alpha = ',sprintf('%.2f',alpha(ii)));
end
 legend(Legenda,'Location','best')
 xlabel('log(N)')
 ylabel('<t_0>/std')
 print(gcf,'mean_std.png','-dpng')


%% creating theoretical curve
h = 0; % external field
modo = 'static'; % random|static
summary_mean = zeros(length(alpha),length(N),1);    
summary_std = summary_mean;

tic;
[summary_mean, summary_std] = generate_histogram(h,alpha,N,modo,time_init,time_end,time_steps);
toc;

figure(1)
show_this(pack_colors, '-', alpha, log(N),log(summary_mean),'mean','log(N)','<t_0>')


figure(2)
show_this(pack_colors, '-', alpha, log(N),log(summary_std),'mean','log(N)','<t_0>')
disp('Task was done')
