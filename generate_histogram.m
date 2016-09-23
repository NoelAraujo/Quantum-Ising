function [summary_mean, summary_std] = generate_histogram(d,J,h,alpha,N,modo,time_init,time_end,time_steps)

%% Simulating

% On these matrices I will store the important data
summary_mean = zeros(length(alpha),length(N),1);    
summary_std = summary_mean;

for aa = 1:length(alpha)
    %clc
    fprintf('Alpha : %d/%d -- %.2f/%.2f \n',aa,length(alpha),alpha(aa),alpha(end))
    for nn=1:length(N)
        %tic;
        fprintf('N : %d/%d -- %d/%d \n',nn,length(N),N(nn),N(end))
        [d_ij, ~, spin0, time_span] = generate_data(modo, N(nn), time_init,time_end,time_steps,alpha(aa),d,J);
        %toc;
        spin = zeros(N(nn),time_steps);
        %tic;
        if strcmp(modo,'static')
            for i=1:1
                logP = 0*time_span;
                for j=1:N(nn)
                    logP = logP + log(abs(cos(2*time_span/(d_ij(i,j).^alpha(aa)))));
                end
                spin(i,:) = spin0(i).*cos(2*h*time_span).*exp(logP);
            end
        else
            for i=1:N(nn)
                logP = 0*time_span;
                for j=1:N(nn)
                    logP = logP + log(abs(cos(2*time_span/(d_ij(i,j).^alpha(aa)))));
                end
                spin(i,:) = spin0(i).*cos(2*h*time_span).*exp(logP);
            end
        end
        %toc;
        %% Fit y = A*exp((-t/t_0)^2) --- nao esquecer da cte
        % --> log(y) = -(t/t_0)^2
        % --> sqrt(-log(y)) = t/t_0
        % --> [[[    t = t_0*sqrt(-log(y))    ]]]
        % t and t_0 are linear as seen at the equation above
        % this equation is used to compute the fit
        %tic;
        t_0 = zeros(1,N(nn));
        if strcmp(modo,'random')
            interval_ii = 1:N(nn);
        else
            interval_ii = 1:1; % on static case, all particle have the same curve
        end
        t_0 = 0;
        for ii=interval_ii
            time_span0 = time_span;
            half_life = abs(spin(ii,:));
            remover = find(half_life < half_life(1)/2);

            time_span0(remover) = [];
            half_life(remover) = [];
            %logP_1 = log(half_life./(0.5.*cos(2*h*time_span0)));
            
            %figure(1)
            %[c, gof] = fit(time_span0'.^2, N(nn)^(2*alpha(aa)-1)*logP_1','poly1');
            %plot(c,time_span0.^2, N(nn)^(2*alpha(aa)-1)*logP_1);
            %drawnow
            %stem(time_span0'.^2,sqrt(-log(half_life))')
            [c, gof] = fit(time_span0',sqrt(-log(half_life/half_life(1)))','poly1');
%             figure(2)
%             plot(c,time_span0'.^1,sqrt(-log(half_life))')
%             drawnow
%             figure(3)
%             plot(time_span0,sqrt(-log(half_life/half_life(1))))
%             pause
            t_0(ii) = 1/c.p1;
            %t_0(ii) = gof.rsquare;
        end
        %toc;
        [summary_mean(aa,nn,1), summary_std(aa,nn,1)] = normfit(t_0); % get the mean and sigma to plot on the histogram
        %summary_mean(aa,nn,1) = mu;
        %summary_std(aa,nn,1) = sigma;
        
    end
end
