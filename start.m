clear all; close all; clc;

%% Initial parameters
d = 1; % dimension
J = 1; 
h = 0; % external field
alpha = [0];
% ref: [Romain Bachelard, Michael Kastner]Universal Threshold for the Dynamical Behavior of Lattice Systems with Long-Range Interactions

N = 100:200:2000;
modo = 'random'; % random|static

time_init = 1e-9;
time_end = 1;
time_steps = 100;

repetitions = 1;

summary_mean = zeros(length(alpha),length(N),1);    
summary_std = summary_mean;

for ii=1:repetitions
    tic;
	fprintf('\n%d/%d-',ii,repetitions)
    [summary_mean_temp, summary_std_temp] = generate_histogram(d,J,h,alpha,N,modo,time_init,time_end,time_steps);
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
for ii=1:length(alpha)
    h = plot(log(N),log(summary_mean(ii,:,1)),strcat(pack_colors(ii),'o-'));
    set(h, 'MarkerFaceColor', get(h, 'Color'));
    hold on
    Legenda{ii} = strcat('\alpha = ',sprintf('%.2f',alpha(ii)));
end
legend(Legenda,'Location','best')
xlabel('N')
ylabel('mean of t_0')
print(gcf,'mean.png','-dpng')

figure(2);
for ii=1:length(alpha)
    h = plot(N,summary_std(ii,:,1),strcat(pack_colors(ii),'o-'));
    set(h, 'MarkerFaceColor', get(h, 'Color'));
    hold on
    Legenda{ii} = strcat('\alpha = ',sprintf('%.2f',alpha(ii)));
end

legend(Legenda,'Location','best')
xlabel('N')
ylabel('std of t_0')
print(gcf,'std.png','-dpng')


figure(3);
for ii=1:length(alpha)
    h = plot(log(N),log(summary_std(ii,:,1)./summary_mean(ii,:,1)),strcat(pack_colors(ii),'o-'));
    set(h, 'MarkerFaceColor', get(h, 'Color'));
    hold on
    Legenda{ii} = strcat('\alpha = ',sprintf('%.2f',alpha(ii)));
end
return
%% creating theoreical curve

h = 0; % external field
modo = 'static'; % random|static
summary_mean = zeros(length(alpha),length(N),1);    
summary_std = summary_mean;

tic;
[summary_mean, summary_std] = generate_histogram(d,J,h,alpha,N,modo,time_init,time_end,time_steps);
toc;

figure(1)
for ii=1:length(alpha)
    h = plot(log(N),log(summary_mean(ii,:,1)),strcat(pack_colors(ii),'o-'));
    hold on
    Legenda{ii} = strcat('\alpha = ',sprintf('%.2f',alpha(ii)));
end
print(gcf,'mean.png','-dpng')



figure(2)
for ii=1:length(alpha)
    h = plot(N,summary_std(ii,:,1),strcat(pack_colors(ii),'o--'),'Linewidth',3);
    hold on
    Legenda{ii} = strcat('\alpha = ',sprintf('%.2f',alpha(ii)));
end
print(gcf,'std.png','-dpng')
disp('Task was done')
