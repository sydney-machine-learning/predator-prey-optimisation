
%%parameters
T =1000;
Np =30;
D = 100;
F_index=14;
[down,up]=test_functions_range(Findex);
lb = down*ones(1,D);
ub = up*ones(1,D);
%%
prob =@test_functions;
f = NaN(Np,1);                      % Vector to store the fitness function value of the population members

BestFitIter = NaN(T+1,1);           % Vector to store the best fitness function value in every iteration

%% Particle Swarm Optimization

P = repmat(lb,Np,1) + repmat((ub-lb),Np,1).*rand(Np,D);   % Generation of the initial population
v = repmat(lb,Np,1) + repmat((ub-lb),Np,1).*rand(Np,D);   % Generation of the initial population

for p = 1:Np
    f(p) = prob(P(p,:),F_index);            % Evaluating the fitness function of the initial population
end

pbest = P;                          % Initialize the personal best solutions
f_pbest = f;                        % Initialize the fitness of the personal best solutions


[f_gbest,ind] = min(f_pbest);       % Determining the best objective function value
gbest = P(ind,:);                   % Determining the best solution

BestFitIter(1) = f_gbest;
c1 =2;c2=2;w=0.99;                  % Hyperparameters
for t = 1:T
             
    for p = 1:Np
        
        v(p,:) = w*v(p,:) + c1*rand(1,D).*(pbest(p,:)-P(p,:)) + c2*rand(1,D).*(gbest - P(p,:)); % generate new velocity
        
        P(p,:) = P(p,:) + v(p,:);   % Update the position (generate new solution)
        
        P(p,:) = max(P(p,:),lb);    % Bounding the violating variables to their lower bound
        P(p,:) = min(P(p,:),ub);    % Bounding the violating variables to their upper bound
        
        f(p) = prob(P(p,:),F_index);        % Determining the fitness of the new solution
        
        if f(p) < f_pbest(p)
            
            f_pbest(p) = f(p);      % updating the fitness function value of the personal best solution
            pbest(p,:) =  P(p,:);   % updating the personal best solution
            
            if f_pbest(p)< f_gbest
                
                f_gbest = f_pbest(p); % updating the fitness function value of the global best solution
                gbest = pbest(p,:);   % updating the global best solution
                
            end
            
        end
    end
    
    BestFitIter(t+1) = f_gbest;      % Storing the best value of each iteration
end

bestfitness = f_gbest;
bestsol = gbest;
plot(0:T,(BestFitIter),'*')          
xlabel('Iteration')
ylabel('Best fitness function value')
