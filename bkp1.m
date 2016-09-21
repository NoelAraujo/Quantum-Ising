function [summary_mean, summary_std] = generate_histogram(d,J,h,alpha,N,modo,time_init,time_end,time_steps)

%% Simulating
% !! do not change the nps value, there is an issue !!
nps = N; %nps = Number of Points to Show

% On these matrices I will store the important data
summary_mean = zeros(length(alpha),length(N),1);    
summary_std = summary_mean;

for aa = 1:length(alpha)
    %clc
    fprintf('Alpha : %d/%d -- %.2f/%.2f \n',aa,length(alpha),alpha(aa),alpha(end))
    for nn=1:length(N)
        [d_ij, mouse_idx,y_t,spin0,time_span] = generate_data(modo, N(nn), nps(nn),time_init,time_end,time_steps,alpha(aa),d,J);
        mouse_points = zeros(time_steps,nps(nn));

        fprintf('Particle : %d/%d -- %.d/%.d \n',nn,length(N),N(nn),N(end))
       % tic;
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
        %toc;

        %% Fit y = exp((-t/t_0)^2)
        % --> log(y) = -(t/t_0)^2
        % --> sqrt(-log(y)) = t/t_0
        % --> [[[    t = t_0*sqrt(-log(y))    ]]]
        % t and t_0 are linear as seen at the equation above
        % this equation is used to compute the fit
        for ii=1:nps
            time_span0 = time_span';
            half_life = abs(mouse_points(:,ii));
            time_span0(find(half_life < 0.25)) = [];
            half_life(find(half_life < 0.25)) = [];
            %[c, gof] = fit(time_span0,(-log(half_life)).^0.5,'poly1');
            %t_0(ii) = c.p1;
            %t_0(ii) = gof.rsquare;
            plot( (time_span0).^2,(-log(half_life)).^1,'o')
            pause(1)
            close all
        end
        
        [mu, sigma] = normfit(t_0); % get the mean and sigma to plot on the histogram
        summary_mean(aa,nn,1) = mu;
        summary_std(aa,nn,1) = sigma;

    end

    
end

end