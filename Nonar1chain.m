% Non-AR(1) Markov Chain
clear all;clc;

%% Part A
epsilon = logninv([0.1, 0.3, 0.5, 0.7, 0.9],1)';
pi = ones(5,1)*0.2;
%h = histfit(epsilon, 5,'lognormal');

%% Part B
% To define a markov transition matrix, note that diagonal elements
% will be \gamma = 0.9 + 0.1 * 0.2 and the non-diagonal will be 0.1*0.2
Pi = ones(5,5)*0.2*0.1;
for i = 1:5
   Pi(i,i) = 0.9+0.1*0.2;
end

% Next, we will create such a matrix in a numerical form. For this,
% we will simulate a path of states using the probabilities given above
N = 99999999;
states = ones(N, 1);
permutation = randperm(5);
states(1) = permutation(1);

tic
% Simulation of path
for j = 2:N
    probability = rand(1);
    if (0 <= probability)&&(probability <= 0.9)
        states(j) = states(j-1);
    else 
        permutation = randperm(5);
        states(j) = permutation(1);
    end
end

% Construct transition matrix from the path
y = zeros(5, 1);
nPi = zeros(5,5);
for k=1:N-1
    y(states(k)) = y(states(k)) + 1;
    nPi(states(k),states(k+1)) = nPi(states(k),states(k+1)) + 1;
end
nPi = bsxfun(@rdivide,nPi,y); 
toc

% 56 seconds in my computer...

%% Part D
% By the law of large numbers, we can think of the initial distribution
% as the proportion of individuals in a given state, and thus just iterate,
% starting at f_1,0 = 1 and f_i,0 = 0, until convergence. 

f = [1, 0, 0, 0, 0]';
tol = 1E-6;
maxit = 150;
it = 1;
distributions = zeros(5, maxit);
distributions(:,1) = f;
dif=1;

while (max(abs(dif)) > tol) && (it ~= maxit)
    distributions(:,it+1) = Pi*distributions(:,it);
    dif = distributions(:,it+1) - distributions(:,it);
    it = it +1;
end


distributions(:,it)

