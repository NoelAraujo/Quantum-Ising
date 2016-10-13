% If I want theoretical curves, I may uncomment this line
%function [summary_mean, summary_std] = generate_histogram(d,J,h,alpha,N,modo,time_init,time_end,time_steps)
function [summary_mean, summary_std] = generate_histogram(h,alpha,N,modo,time_init,time_end,time_steps)

% The goals of this file are :
% - Generate a lot of spins profiles
% - Get its decaying constant
% - Store them on structured order
% - Return only the proced info to the user


% On these matrices I will store the important data
summary_mean = zeros(length(alpha),length(N),1);
summary_std = summary_mean;

for aa = 1:length(alpha)
    %fprintf('Alpha : %d/%d -- %.2f/%.2f \n',aa,length(alpha),alpha(aa),alpha(end))
    for nn=1:length(N)
        
        %fprintf('N : %d/%d -- %d/%d \n',nn,length(N),N(nn),N(end))
        
        % If I need to theoretical curve, I must to change the
        % generate_data-function. All that I have to do is to comments some
        % parts and uncomment others
        %[d_ij,y_t,spin0,time_span] = generate_data(modo, N(nn),time_init,time_end,time_steps,alpha(aa),d,J);
        
        [d_ij,spin0,time_span] = generate_data(modo, N(nn), time_init,time_end,time_steps);
        spin = zeros(N(nn),time_steps);
        
        % the following loop is written only to be crystal clear to understand, it
        % can be improved
        if strcmp(modo,'static')
            for i=1:1 %only 1 loop because all particle have the same curve
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
            % plot(time_span,teorico,'k','LineWidth',3); hold on;
            % plot(time_span,abs(spin(i,:)) ); drawnow
        end
        
        %% Fit y = y0*exp((-t/t_0)^2) --- DO NOT FORGET THE CONSTANT
        % --> log(y/y0) = -(t/t_0)^2
        % --> sqrt(-log(y/y0)) = t/t_0
        % --> [[[    t = t_0*sqrt(-log(y/y0))    ]]]
        % t and t_0 are linear as seen at the equation above
        % this equation is used to compute the fit
        
        t_0 = zeros(1,N(nn));
        if strcmp(modo,'random')
            interval_ii = 1:N(nn);
        else
            interval_ii = 1:1; % on static case, all particle have the same curve
        end
        
        parfor ii=interval_ii
            time_span0 = time_span;
            half_life = abs(spin(ii,:));
            % My goal is to get the time to decaying HALF of the spin, then
            % I remove all the points below this level to perform a fit and
            % get is time to decay
            remover = find(half_life < half_life(1)/2);
            time_span0(remover) = [];
            half_life(remover) = [];
            
            [c, gof] = fit(time_span0',sqrt(-log(half_life/half_life(1)))','poly1');
            % There are some decaying behavior that looks like damped
            % oscilations. Then the linearization that I did on the
            % previous line is not good. Then I apply a specialized fit.
            % Look on the validation files for more information
            % --> There are some behaviors that also occur, but they are
            % harder to deal, then I will simply take then on the
            % calculation, and hope that they not f#ck my final result.<--
            if gof.rsquare < 0.9
                t_0(ii) = 1/damped_oscilations_coeff(time_span, spin(ii,:));
            else
                t_0(ii) = 1/c.p1;    
                %t_0(ii) = gof.rsquare; % For debugging
            end
        end
        
        % get the mean and sigma that I would see IF I plot the histogram of the data
        [summary_mean(aa,nn,1), summary_std(aa,nn,1)] = normfit(t_0); 
        
    end
end
