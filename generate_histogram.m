function [summary_mean, summary_std] = generate_histogram(d,J,h,alpha,N,modo,time_init,time_end,time_steps)

%% Simulating

% On these matrices I will store the important data
summary_mean = zeros(length(alpha),length(N),1);    
summary_std = summary_mean;
spin = zeros(N,time_steps);
for aa = 1:length(alpha)
    %clc
    fprintf('Alpha : %d/%d -- %.2f/%.2f \n',aa,length(alpha),alpha(aa),alpha(end))
    for nn=1:length(N)
        %tic;
        [d_ij, y_t, spin0, time_span] = generate_data(modo, N(nn), time_init,time_end,time_steps,alpha(aa),d,J);
        %toc;
        %tic;
        for i=1:N
            logP = 0*time_span;
            for j=1:N
                logP = logP + log(abs(cos(2*time_span/(d_ij(i,j).^alpha(aa)))));
            end
            spin(i,:) = spin0(i).*cos(2*h*time_span).*exp(logP);
        end
        %toc;
        %% Fit y = exp((-t/t_0)^2)
        % --> log(y) = -(t/t_0)^2
        % --> sqrt(-log(y)) = t/t_0
        % --> [[[    t = t_0*sqrt(-log(y))    ]]]
        % t and t_0 are linear as seen at the equation above
        % this equation is used to compute the fit
        %tic;
        parfor ii=1:N
            time_span0 = time_span;
            half_life = abs(spin(ii,:));
            remover = find(half_life < 0.10);

            time_span0(remover) = [];
            half_life(remover) = [];
            %logP_1 = log(half_life./(0.5.*cos(2*h*time_span0)));
            %plot(time_span0.^2, N(nn)^(2*alpha(aa)-1)*logP_1);
            %drawnow
            %[c, gof] = fit(time_span0',logP_1','poly1');
            [c, gof] = fit(time_span0',half_life','poly1');
            t_0(ii) = c.p1;
            t_0(ii) = gof.rsquare;
        end
        %toc;
        [mu, sigma] = normfit(t_0); % get the mean and sigma to plot on the histogram
        summary_mean(aa,nn,1) = mu;
        summary_std(aa,nn,1) = sigma;

    end
end