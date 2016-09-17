clear all; close all; clc;

%% Initial parameters
d = 1; % dimension
J = 1; 
h = 0; % external field
alpha = [0.2 0.4 0.6 0.8];
% ref: [Romain Bachelard, Michael Kastner]Universal Threshold for the Dynamical Behavior of Lattice Systems with Long-Range Interactions

N = [20 40 60 80 100];
modo = 'static'; % random|static

time_init = 1e-8;
time_end = 1;
time_steps = 100;

repetitions = 10;

summary_mean = zeros(length(alpha),length(N),1);    
summary_std = summary_mean;

for ii=1:repetitions
    tic;
	fprintf('%d/%d-',ii,repetitions)
    [summary_mean_temp, summary_std_temp] = generate_histogram(d,J,h,alpha,N,modo,time_init,time_end,time_steps);
    summary_mean = summary_mean + summary_mean_temp;
    summary_std = summary_std + summary_std_temp;
    toc;
end

summary_mean = summary_mean./repetitions;
summary_std = summary_std./repetitions;
%%
close all
figure;
for ii=1:length(alpha)
    h = plot(N,summary_mean(ii,:,1),'o-');
    set(h, 'MarkerFaceColor', get(h, 'Color'));
    hold on
    Legenda{ii} = strcat('\alpha = ',sprintf('%.2f',alpha(ii)));
end
legend(Legenda,'Location','best')
xlabel('N')
ylabel('mean of t_0')
print(gcf,'mean.png','-dpng')

figure;
for ii=1:length(alpha)
    h = plot(N,summary_std(ii,:,1),'o-');
    set(h, 'MarkerFaceColor', get(h, 'Color'));
    hold on
    Legenda{ii} = strcat('\alpha = ',sprintf('%.2f',alpha(ii)));
end

legend(Legenda,'Location','best')
xlabel('N')
ylabel('std of t_0')
print(gcf,'std.png','-dpng')

disp('Task was done')