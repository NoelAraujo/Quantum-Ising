clear all; close all; clc;

%% Initial parameters
d = 1; % dimension
J = 1; 
h = 0; % external field
alpha = 0.65;
% ref: [Romain Bachelard, Michael Kastner]Universal Threshold for the Dynamical Behavior of Lattice Systems with Long-Range Interactions

N = 250;
modo = 'random'; % random|static
nps = 250; %nps = Number of Points to Show

time_init = 0;
time_end = 1.5;
time_steps = 100;



%% avoid edit below this line

if strcmp('random',modo)
    % Create N points over a line with random distance
    distance = 1;
    pd = makedist('Uniform','lower',0,'upper',(N*distance));
    position = random(pd,[1 N]);
    %figure(1);histogram(position)
elseif strcmp('static',modo)
    % Create N points over a line with fixed unit distance
    distance = 1;
    position = distance:distance:(N*distance);
    % all the results should be equal fot the static case
else
    disp('!! Problems !!')
    disp('The end')
    return
end

if nps>N
    disp('!! Problems with nps !!')
    return
end


% the mouse indices are 'nps' random points that I follow their behavior along the time
    % I'm using the word "mouse" in reference to the mouse at biologic
    % experiments done with mouses
    if nps == N
        mouse_idx = 1:N;
    else
        mouse_idx = randi([1 N],[1 nps]); 
    end

%% Create a random 1/2-spin for each particle
spin0 = 0.5*randsrc(1,N);

%% Formula 2
% get the correct distance between all particles
r = position; % I've created the variable "r" to avoid change notation. The following code I've used on other projects
d_ij1 = sqrt(bsxfun(@plus,dot(r,r,1)',dot(r,r,1))-2*(r'*r)); %position of each particle with each particle
d_ij2 = (N*distance) - d_ij1; % total lenght - distance between particles

d_ij = min(d_ij1,d_ij2); % In periodic boundary condition, I need the smallest distance.

d_ij(find(d_ij==0)) = -10;% trick to speed up the program
%     mouse_idx = 1:20;
% load d_ij.mat
%% Black Curve
mouse_points = zeros(time_steps,nps);
y_t = zeros(1,time_steps);
time_span = linspace(time_init,time_end,time_steps);

if alpha >= d/2
    for t = 1:time_steps
        y_t(t) = 0.5*exp(-(2^(1+d-2*alpha)*(pi)^(d/2) )/((2*alpha-d)*gamma(d/2))*(4*J*time_span(t)/pi)^(d/alpha));
    end
else
    for t = 1:time_steps
        y_t(t) = 0.5*exp(-(2^(5+2*alpha-d)*(pi)^(d/2-2) )/((d-2*alpha)*gamma(d/2))*(J^2*time_span(t)^2)*(N)^(1-(2*alpha)/d));
    end
end

%% Simulation
tic;
for t = 1:time_steps
   % disp(time_span(t))
    spin = spin0; 
        
    for i=1:N
        spin_temp = spin(i);% trick to speed up simulation
        
        for j=1:N
                spin_temp = spin_temp*cos((2*J*time_span(t))./abs(d_ij(i,j).^alpha));
        end
        spin(i) = spin_temp/cos( (2*J*time_span(t))/abs(-10^alpha) ); % trick to speed up simulation
    end
    
    spin = spin*cos(2*h*time_span(t)); % on the for-loops I've compute the multiplications of all cos(). Now I multiply by the initial spin value
    mouse_points(t,:) = spin(mouse_idx);
    
end
toc;

%% Plot only some points
for ii=1:nps
    plot(time_span, (abs(mouse_points(:,ii))),'color',rand([1 3]))
    hold on
end
    plot(time_span, (y_t),'-k','LineWidth',5)