clear all; close all; clc;

%% Initial parameters
d = 1; % dimension
J = 1; 
h = 0; % external field
alpha = [0.1 0.15 0.6 0.9];
% ref: [Romain Bachelard, Michael Kastner]Universal Threshold for the Dynamical Behavior of Lattice Systems with Long-Range Interactions

N = [20 50 100 150];
modo = 'random'; % random|static

time_init = 0;
time_end = 1;
time_steps = 100;




%% Simulating
%     system('rd /s /q C:\Users\Simus\Desktop\Noel\Quantum-Ising-master\Results\alpha')
%     system('mkdir C:\Users\Simus\Desktop\Noel\Quantum-Ising-master\Results\alpha')
%     system('del C:\Users\Simus\Desktop\Noel\Quantum-Ising-master\*.png')

% system('rm -r Results/alpha')
% system('mkdir Results/alpha')
% system('rm *.png')

% !! do not change the nps value, there is an issue !!
nps = N; %nps = Number of Points to Show

% On these matrices I will store the important data
summary_mean = zeros(length(alpha),length(N),1);    
summary_std = summary_mean;

for aa = 1:length(alpha)
    clc
    fprintf('Alpha : %d/%d -- %.2f/%.2f \n',aa,length(alpha),alpha(aa),alpha(end))
    %pause(0.5)
    for nn=1:length(N)
        [d_ij, mouse_idx,y_t,spin0,time_span] = generate_data(modo, N(nn), nps(nn),time_init,time_end,time_steps,alpha(aa),d,J);
        mouse_points = zeros(time_steps,nps(nn));

        fprintf('Particle : %d/%d -- %.d/%.d \n',nn,length(N),N(nn),N(end))
        tic;
        for t = 1:time_steps
           % disp(time_span(t))
            spin = spin0; 
            cos_temp = cos((2*J*time_span(t))./abs(d_ij.^alpha(aa)));% trick to speed up simulation
            for i = 1:N(nn)
                spin_temp = spin(i);% trick to speed up simulation
                for j=1:N(nn)
                         spin_temp = spin_temp*cos_temp(i,j);
                end
                spin(i) = spin_temp/cos( (2*J*time_span(t))/abs(-10^alpha(aa)) ); % trick to speed up simulation
            end
            spin = spin*cos(2*h*time_span(t)); % on the for-loops I've compute the multiplications of all cos(). Now I multiply by the initial spin value
            mouse_points(t,:) = spin(mouse_idx);
        end
        toc;

        %% Plot only some points
%         fig = figure;
%         for ii=1:nps
%             plot(time_span, (abs(mouse_points(:,ii))),'color',rand([1 3]))
%             hold on
%         end
%         plot(time_span, (y_t),'-k','LineWidth',5)
%         xlabel('Time')
%         ylabel('< Spin >')
%         titulo = sprintf('alpha=%.2f, N=%d',alpha(aa),N(nn));
%         title(titulo)

%         set(gcf, 'renderer', 'opengl')
%         print(fig,strcat(titulo,'.png'),'-dpng')
%         close all
%         
        %% Fit y = exp((-t/t_0)^2)
        % --> sqrt(-log(y)) = t/t_0
        % --> [[[    t = t_0*sqrt(-log(y))    ]]]
        % t and t_0 are linear as seen at the equation above
        % this equation is used to compute the fit
        for ii=1:nps
            time_span0 = time_span';
            half_life = abs(mouse_points(:,ii));
            time_span0(find(half_life < 0.25)) = [];
            half_life(find(half_life < 0.25)) = [];
            [c, gof] = fit(time_span0,sqrt(-log(half_life)),'poly1');
            t_0(ii) = c.p1;
        end
        
        [mu, sigma] = normfit(t_0); % get the mean and sigma to plot on the histogram
        summary_mean(aa,nn,1) = mu;
        summary_std(aa,nn,1) = sigma;
%         titulo = sprintf('alpha = %.2f,N = %d, mean = %.2f, sigma = %.3f',alpha(aa),N(nn),mu,sigma);
%         histfit(t_0); %create histograma with norm fit automatically
%         title(titulo)
%         print(gcf,strcat(titulo,'.png'),'-dpng')
%         close all
    end
    %pause(1)

%     system(strcat('mkdir Results/alpha/',sprintf('%.2f',alpha(aa))))
%     system(strcat('mv *png Results/alpha/',sprintf('%.2f',alpha(aa))))

%     system(strcat('mkdir C:\Users\Simus\Desktop\Noel\Quantum-Ising-master\Results\alpha\',sprintf('%.2f',alpha(aa))))
%     system(strcat('move C:\Users\Simus\Desktop\Noel\Quantum-Ising-master\*.png C:\Users\Simus\Desktop\Noel\Quantum-Ising-master\Results\alpha\',sprintf('%.2f\\',alpha(aa))))
    
end

figure;
ax = axes;
imagesc(alpha,N,summary_mean,'Parent',ax);
colorbar
title('Mean(alpha,N)')
xlabel('alpha')
ylabel('N')
set(ax,'YDir','normal')
xlim([min(alpha) max(alpha)])
ylim([min(N) max(N)])
print(gcf,'mean.png','-dpng')

figure;
ax = axes;
imagesc(alpha,N,summary_std,'Parent',ax);
colorbar
title('STD(alpha,N)')
xlabel('alpha')
ylabel('N')
set(ax,'YDir','normal')
xlim([min(alpha) max(alpha)])
ylim([min(N) max(N)])
print(gcf,'std.png','-dpng')

disp('Task was done')