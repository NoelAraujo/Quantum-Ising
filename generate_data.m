%function [d_ij,y_t,spin0,time_span] = generate_data(modo, N, time_init,time_end,time_steps,alpha,d,J)
function [d_ij,spin0,time_span] = generate_data(modo, N, time_init,time_end,time_steps) 
% The goals of this file are:
% - Create N random points with one dimensional position (1D)
% - There position may be random or uniformily distributed
% - Compute the mean spin of the system
% - Create theoretical curve if needed

    if strcmp('random',modo)
        % Create N points over a line with random distance
        distance = 1;
        pd = makedist('Uniform','lower',0,'upper',(N*distance));
        position = random(pd,[1 N]);
        %figure(1);histogram(position)
    elseif strcmp('static',modo)
        % Create N points over a line with fixed unit distance
        distance = 1;
        position = 0:distance:((N-1)*distance);
        % all the results should be equal fot the static case
    else
        disp('!! Problems !!')
        disp('The end')
        return
    end


    %% Create a random [-1,+1]-spin for each particle
    spin0 = randsrc(1,N);

    %% Formula 2
    % get the correct distance between all particles
    r = position; % I've created the variable "r" to avoid change notation. The following code I've used on other projects
    d_ij1 = sqrt(bsxfun(@plus,dot(r,r,1)',dot(r,r,1))-2*(r'*r)); %position of each particle with each particle
    
    % Trying something new here
%     to_be_used = round(sqrt(N));
%     out = 0;
% 
%     while out == 0
%         index_to_be_used = randi([1 N],[1 to_be_used]);
%         if sum(diff(sort(index_to_be_used))==0) == 0
%            out = 1;
%         end
%     end
%     d_ij1 = zeros(to_be_used,N);
%     for k = 1:to_be_used
%         d_ij1(k,:) = abs(position(index_to_be_used(k)) - position);
%     end
    d_ij2 = (N*distance) - d_ij1; % total lenght - distance between particles

    d_ij = min(d_ij1,d_ij2); % In periodic boundary condition, I need the smallest distance.
    d_ij(find(d_ij==0)) = 0.5;% trick to speed up the program
        
    %% This is necessary to overcome the problem of very fast oscillations 
    % They decay so fast the make problems during the simulation
    % This make more points near the origin, then I can see its behavior
    % more precisily before it goes to zero - thus I can estimate is half
    % life time
    
    %time_span = linspace(time_init,time_end,time_steps); % Original linspace
    
    time_span = logspace(log10(time_init),log10(time_end),time_steps);
    %% When I need the theoretical cuve I may use this section
    % If not in use, this section should be in comments to performance
%     y_t = zeros(1,time_steps);
%     if alpha >= d/2
%         for t = 1:time_steps
%             y_t(t) = exp(-(2^(1+d-2*alpha)*(pi)^(d/2) )/((2*alpha-d)*gamma(d/2))*(4*J*time_span(t)/pi)^(d/alpha));
%         end
%     else
%         for t = 1:time_steps
%             y_t(t) = exp(-(2^(5+2*alpha-d)*(pi)^(d/2-2) )/((d-2*alpha)*gamma(d/2))*(J^2*time_span(t)^2)*(N)^(1-(2*alpha)/d));
%         end
%     end


end
