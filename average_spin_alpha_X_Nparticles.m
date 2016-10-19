clear all; close all; clc;
% Goals of this file are:
%   - Vary alpha and N
%   - Generate plots of average time(<t0>) for decaying spins
%   - On a plot log-log, there is a linear relationship between <t0> and
%   the number of particles
%% Initial parameters
d = 1; % dimension
J = 1; 
h = 0; % external field
alpha = [0.2 0.4 0.6 0.8];
% ref: [Romain Bachelard, Michael Kastner]Universal Threshold for the Dynamical Behavior of Lattice Systems with Long-Range Interactions

N = round(linspace(100,500,20));
modo = 'static'; % random|static

time_init = 1e-9;
time_end = 1;
time_steps = 100;

repetitions = 3;
 
summary_mean = zeros(length(alpha),length(N),1);    
summary_std = summary_mean;

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



%% Plots
pack_colors = ['r','b','k','g'];


close all
figure(1);
show_this(pack_colors, 'o-', alpha, log(N),log(summary_mean),'mean','log(N)','log(<t_0>)')

figure(2);
show_this(pack_colors, 'o-',alpha, log(N),log(summary_std),'sigma','log(N)','log(std)')

%figure(3);
% for ii=1:length(alpha)
%     h = plot(log(N),log(summary_std(ii,:,1)./summary_mean(ii,:,1)), strcat(pack_colors(ii),'o-'));
%     set(h, 'MarkerFaceColor', get(h, 'Color'));
%     hold on
%     Legenda{ii} = strcat('\alpha = ',sprintf('%.2f',alpha(ii)));
% end
%  legend(Legenda,'Location','best')
%  xlabel('log(N)')
%  ylabel('<t_0>/std')
% print(gcf,'mean_std.png','-dpng')


%% creating theoretical curve
h = 0; % external field
modo = 'static'; % random|static
summary_mean = zeros(length(alpha),length(N),1);    
%summary_std = summary_mean;

tic;
[summary_mean, ~] = generate_histogram(h,alpha,N,modo,time_init,time_end,time_steps);
toc;

figure(1)
show_this(pack_colors, '-', alpha, log(N),log(summary_mean),'mean','log(N)','<t_0>')


%figure(2)
%show_this(pack_colors, '-', alpha, log(N),log(summary_std),'mean','log(N)','<t_0>')


disp('Task was done')
