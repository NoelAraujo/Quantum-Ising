clear all; close all; clc;

%% Initial parameters
d = 1; % dimension
J = 1; 
h = 0; % external field
alpha = 0.1:0.1:0.9;
% ref: [Romain Bachelard, Michael Kastner]Universal Threshold for the Dynamical Behavior of Lattice Systems with Long-Range Interactions

N = 250:250:3000;
modo = 'random'; % random|static
nps = N; %nps = Number of Points to Show

time_init = 0;
time_end = 1;
time_steps = 100;




%% Simulating
for aa = 1:length(alpha)
    for nn=1:length(N)
            [d_ij, mouse_idx,y_t,spin0,time_span] = generate_data(modo, N(nn), nps(nn),time_init,time_end,time_steps,alpha(aa),d,J);
            mouse_points = zeros(time_steps,nps(nn));

            disp('Working on')
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
        fig = figure;
        for ii=1:nps
            plot(time_span, (abs(mouse_points(:,ii))),'color',rand([1 3]))
            hold on
        end
            plot(time_span, (y_t),'-k','LineWidth',5)
            xlabel('Time')
            ylabel('< Spin >')
            titulo = sprintf('alpha=%.2d, N=%d',alpha(aa),N(nn));
            title(titulo)

            set(gcf, 'renderer', 'opengl')
            print(fig,strcat(titulo,'.png'),'-dpng')
            close all
    end
    system(strcat('rm -r Results/',sprintf('alpha=%.2d',alpha(aa))))
    system(strcat('mkdir Results/',sprintf('alpha=%.2d',alpha(aa))))
    system(strcat('mv *png Results/',sprintf('alpha=%.2d',alpha(aa))))
end
