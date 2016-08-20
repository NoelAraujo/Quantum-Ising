clear all; close all; clc;

%% Create N points over a line with fixed distance
N = 500;
distance = 1;
position = distance:distance:(N*distance);

%% Create a random 1/2-spin for each particle
spin = 0.5*randsrc(1,N);

%% Initial parameters ( formaula 2 )
d = 1; % dimension
J = 1; 
h = 0;
alpha = 1;
% ref: [Romain Bachelard, Michael Kastner]Universal Threshold for the Dynamical Behavior of Lattice Systems with Long-Range Interactions

%% Formula (2)
% get the correct distance between all particles
r = position;
d_ij1 = sqrt(bsxfun(@plus,dot(r,r,1)',dot(r,r,1))-2*(r'*r)); %position of each particle with each particle
d_ij2 = (N*distance) - d_ij1; % total lenght - distance between particles

d_ij = min(d_ij1,d_ij2); % In periodic boundary condition, I need the smallest distance.


for t = 1
    for i=1:N
        for j=1:N
            if i~=j
                spin(i) = spin(i)*cos(2*h*t)*cos( (2*J*t)/(d_ij(i,j)^alpha) );
            end
        end
    end
end
stem(spin)