clear all; close all; clc;

N = 1000;
modo = 'random'; % random|static


%% avoid edit below this line
if strcmp('random',modo)
    % Create N points over a line with random distance
    distance = 1;
    pd = makedist('Uniform','lower',0,'upper',(N*distance));
    position = random(pd,[1 N]);
    histogram(position)
elseif strcmp('static',modo)
    % Create N points over a line with fixed unit distance
    distance = 1;
    position = distance:distance:(N*distance);    
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
alpha = 0.45;
% ref: [Romain Bachelard, Michael Kastner]Universal Threshold for the Dynamical Behavior of Lattice Systems with Long-Range Interactions

%% Formula 2
% get the correct distance between all particles
r = position; % I've created the variable "r" to avoid change notation. The following code I've used on other projects
d_ij1 = sqrt(bsxfun(@plus,dot(r,r,1)',dot(r,r,1))-2*(r'*r)); %position of each particle with each particle
d_ij2 = (N*distance) - d_ij1; % total lenght - distance between particles

d_ij = min(d_ij1,d_ij2); % In periodic boundary condition, I need the smallest distance.

%% Simulation
for t = linspace(0.01,1,100)
    disp(t)
    spin = spin0; 
    
    if alpha >= d/2
        y=0*spin+0.5*exp(-(2^(1+d-2*alpha)*(pi)^(d/2) )/((2*alpha-d)*gamma(d/2))*(4*J*t/pi)^(d/alpha));
    else
        y=0*spin+0.5*exp(-(2^(5+2*alpha-d)*(pi)^(d/2-2) )/((d-2*alpha)*gamma(d/2))*(J^2*t^2)*(N)^(1-(2*alpha)/d));
    end
    
    for i=1:N
        for j=1:N
            if i~=j
                spin(i) = spin(i)*cos( (2*J*t)/abs(d_ij(i,j)^alpha) );
            end
        end
    end
    spin = spin*cos(2*h*t); % on the for-loops I've compute the multiplications of all cos(). Now I multiply by the initial spin value
    stem(spin,'filled','LineStyle','none')
    hold on
    plot(position,y,'r')
    plot(position,-y,'r')
    hold off
    drawnow
end
