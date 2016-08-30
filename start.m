clear all; close all; clc;

N = 1500;
modo = 'random'; % random|static


%% avoid edit below this line
if strcmp('random',modo)
    % Create N points over a line with random distance
    distance = 1;
    pd = makedist('Uniform','lower',0,'upper',(N*distance));
    position = random(pd,[1 N]);
    figure(1);histogram(position)
    % the mouse indices are 10 random points that I follow their behavior along the time
    % I'm using the word "mouse" in reference to the mouse at biologic
    % experiments done with mouses
    mouse_idx = randi([1 N],[1 10]); 
    
elseif strcmp('static',modo)
    % Create N points over a line with fixed unit distance
    distance = 1;
    position = distance:distance:(N*distance);
    % all the results should be equal fot the static case
    mouse_idx = randi([1 N],[1 10]); 
else
    disp('!! Problems !!')
    disp('The end')
    return
end
%% Create a random 1/2-spin for each particle
spin0 = 0.5*randsrc(1,N);

%% Initial parameters ( formaula 2 )
d = 1; % dimension
J = 1; 
h = 0; % external field
alpha = 0.35;
% ref: [Romain Bachelard, Michael Kastner]Universal Threshold for the Dynamical Behavior of Lattice Systems with Long-Range Interactions

%% Formula 2
% get the correct distance between all particles
r = position; % I've created the variable "r" to avoid change notation. The following code I've used on other projects
d_ij1 = sqrt(bsxfun(@plus,dot(r,r,1)',dot(r,r,1))-2*(r'*r)); %position of each particle with each particle
d_ij2 = (N*distance) - d_ij1; % total lenght - distance between particles

d_ij = min(d_ij1,d_ij2); % In periodic boundary condition, I need the smallest distance.

%% Simulation
% figure(2)
mouse_points = zeros(100,10);
y_t = zeros(1,100);
time_counter = 1;
for t = linspace(0.01,0.5,100)
    disp(t)
    spin = spin0; 
    
    if alpha >= d/2
        y=0*spin+0.5*exp(-(2^(1+d-2*alpha)*(pi)^(d/2) )/((2*alpha-d)*gamma(d/2))*(4*J*t/pi)^(d/alpha));
        y_t(time_counter) = 0.5*exp(-(2^(1+d-2*alpha)*(pi)^(d/2) )/((2*alpha-d)*gamma(d/2))*(4*J*t/pi)^(d/alpha));
    else
        y=0*spin+0.5*exp(-(2^(5+2*alpha-d)*(pi)^(d/2-2) )/((d-2*alpha)*gamma(d/2))*(J^2*t^2)*(N)^(1-(2*alpha)/d));
        y_t(time_counter) = 0.5*exp(-(2^(5+2*alpha-d)*(pi)^(d/2-2) )/((d-2*alpha)*gamma(d/2))*(J^2*t^2)*(N)^(1-(2*alpha)/d));
    end
    
    for i=1:N
        for j=1:N
            if i~=j
                spin(i) = spin(i)*cos( (2*J*t)/abs(d_ij(i,j)^alpha) );
            end
        end
    end
    spin = spin*cos(2*h*t); % on the for-loops I've compute the multiplications of all cos(). Now I multiply by the initial spin value
    mouse_points(time_counter,:) = spin(mouse_idx);
    time_counter = time_counter+1;
%     figure(2);stem(spin,'filled','LineStyle','none')
%     hold on
%     figure(2);plot(position,y,'r')
%     figure(2);plot(position,-y,'r')
%     hold off
%     drawnow
end

%%
for ii=1:10
    scatter(linspace(0.01,0.5,100), abs(mouse_points(:,ii)),'filled')
    hold on
end
    plot(linspace(0.01,0.5,100), y_t,'-k','LineWidth',5)