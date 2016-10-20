clear all; close all; clc;

%% Initial parameters
d = 1; % dimension
J = 1; 
h = 0; % external field
alpha = linspace(0.1,0.95,10);
% ref: [Romain Bachelard, Michael Kastner]Universal Threshold for the Dynamical Behavior of Lattice Systems with Long-Range Interactions

Nm(1,:) = round(linspace(1e2,5e2,10));
Nm(2,:) = round(linspace(1e3,1.5e3,10));
% %Nm(3,:) = round(linspace(1e4,1e5,10));
% Nm(3,:) = round(linspace(1e5,1e6,10));
modo = 'random'; % random|static

time_init = 1e-9;
time_end = 1;
time_steps = 100;

repetitions = 3;
pack_colors = ['r','b','g'];

for nn=1:2
    N = Nm(nn,:);
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
    q_mean = 0;
    for ii=1:length(alpha)  
        [coeff_mean,~] = fit(log(N'),log(summary_mean(ii,:)'),'poly1');
        %plot(log(N'),log(summary_mean(ii,:)')); hold on
        %[coeff_std,~]  = fit(log(N),log(summary_std(ii,:)));

        q_mean(ii) = coeff_mean.p1;
        %q_std(ii) = coeff_std.p1;

    end

    %% Plots

    figure(100);
    hold on
    %plot(alpha,q_mean,'-o')
    hh = plot(alpha,q_mean,strcat(pack_colors(nn),'o-'));
    set(hh, 'MarkerFaceColor', get(hh, 'Color'));
    hold on
    Legenda{nn} = strcat('N : ',sprintf('[%.d',min(Nm(nn,:))),',',sprintf('%.d',max(Nm(nn,:))),']');
end
legend(Legenda,'Location','best')
%%
figure(100)
hold on

q_theoretical = min(0,[0:0.01:1]-0.5);
plot([0:0.01:1],q_theoretical,'k--','Linewidth',2)
line([0.5 0.5],[-0.55 0])
ylim([-0.55 0.25])