function [bestsol, bestfitness, BestFitIter, P, f] = QCMBO_2( prob, lb, ub, T, Findex)

%% Algorithm parameters
D = length(lb);                                                            %Number of continuous decision variables 
%T = 1000;                                                                 %Iterations
%prob = @Rastrigin;                                                        %Function handle
Popsizes =[30];                                                     %Multiple Populations sizes vector for comparitive evaluation in terms of population width
for i = 1:length(Popsizes)                                                 %length(Popsizes)
Np = Popsizes(i);                                                          %Np = Popsizes(i);
%lb = -1000*ones(1,D);                                                     %Lower Bound on Decision Variables
%ub = +1000*ones(1,D);                                                     %Upper Bound on Decision Variables

Nm = Np/2;                                                                 %Population of Mice
Nc = Np/2;                                                                 %Population of Cats
N  = Np;                                                                   %Total population of Mice and Cats
%% Stating of CMBO
f = NaN(N,1);                                                              %Vector to store fitness function values                                                        
P = repmat(lb,N,1) + repmat((ub-lb),N,1).*rand(N,D);                       %Generation of the initial population
for p = 1:N
    f(p) = prob(P(p,:),Findex);                                                   %Calculating fitness of the population
end
%% Initialising Vectors for storing best sol. and fit. after every iteration
%IterFit = NaN(T,1);
sol = NaN(T,D);
BestFitIter = NaN(T,1);
%% Sorting initial population
[fsorted ,I] = sort(f);                                                    %Sorting the fitness column vector
Psorted =  P(I,:);                                                         %Sorting Population vector as per sorted fitness
%% Iteration loop
for t = 1: T
Set = Nm;
%% Sorting Population based on fitness values
[fsorted ,I] = sort(fsorted);                                              %Sorting the fitness column vector after every iterations
Psorted =  Psorted(I,:);                                                   %Sorting of population based on sorted fitness after every iterations
%% Selection of Population
Pmice = Psorted(1:Nm ,:);                                                  %Mice population
fmice = fsorted(1:Nm ,:);                                                  %fitness of mice
Pcat  = Psorted(Nm+1:N ,:);                                                %Cats population
fcat  = fsorted(Nm+1:N ,:);                                                %fitness of cats

%% Cat Phase
for j = 1:Nc  
   % b=0.99;
    %u  = 0.01+b*rand;
   
    k = randi([1 Set],1,1);                                                %Selecting Mouse to be attacked
    if(Set <= Nm)
       Set = Set -1;                                                       %Decreasing Scope of Mice as cat attacks increases
    end
   % L = randi([1 2],1,1);                                                 %Generating either 1 or 2 randomly for teaching factor
    
    %Updating Cat after attacking on the Mouse
    c1 = rand;
    c2 = rand;                                                      
    A  = (c1*Pmice(k,:) + c2*Pcat(j,:))/(c1+c2);
   
    %Possible Positions of a quantum Cats (*particles)
    NewCat1 = A - (Psorted(k,:)-c2*Pcat(j,:))*log(1/c1);
    NewCat2 = A + (Psorted(k,:)-c2*Pcat(j,:))*log(1/c1);

    %Bounding New Cats Solution
    NewCat1 = min(ub,NewCat1);
    NewCat1 = max(lb,NewCat1);
    NewCat2 = min(ub,NewCat2);
    NewCat2 = max(lb,NewCat2);

    MIN = NewCat2;

    %Greedy Selection
    if(prob(NewCat1, Findex)<prob(NewCat2, Findex))
        MIN = NewCat1;
    end
    %Selecting best position of the particles
    if (prob(MIN, Findex)<fcat(j))
        fsorted(Nm+j) = prob(MIN, Findex);
        Psorted(Nm+j,:) = MIN;
    end
end
%% Mouse Phase
for j = 1:Nm
    %Increasing Haven scope as mice population decreases
    [~,k] = min(fsorted);
  
    Haven = Psorted(k,:);
    %Implementing Exploring factor E
   
    E = randi([-2 2],1,1);
   
    %Updating Mouse after escaping to the Haven
    Newmice = Pmice(j,:) + E*rand(1, D).*( Haven-Pmice(j,:));
    %Bounding New Mouse Solution
    Newmice = min(ub,Newmice);
    Newmice = max(lb,Newmice);
    %Greedy Selection
    if (prob(Newmice, Findex)<fmice(j))
        fsorted(j) = prob(Newmice, Findex);
        Psorted(j,:) = Newmice;
    end
end
%Storing Best Results
[bestfitness ,ind]   = min(fsorted);                                        
BestFitIter(t,:) = bestfitness;                                            %Best fitness of current iteration is stored in
bestsol = Psorted(ind,:);                                                  
sol(t,:) = bestsol;                                                        %Best Solution of current iteration is stored in
end
%% Plotting best fitness vs iterations
hold on
plot(1:T,(BestFitIter(1:T)),'-', 'LineWidth', 6)                                %Using Semilog Plot to get More Detailed variation
xlabel('Iteration')
ylabel('Best fitness function value')

hold on
end
legend('Np=30','Np=50','Np=100');
title('Rosenbrock')
end