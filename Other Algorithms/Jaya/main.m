%% JAYA algorithms
clc
clear all
close all

%% Problem Definition
F_index=14;
pop = 30;               % Population size
var = 100;                 % Number of design variables
maxGen = 1000;            % Maximum number of iterations
[up,down]=test_functions_range(F_index);
mini = down*ones(1,var);  % Lower Bound of Variables
maxi = up*ones(1,var);   % Upper Bound of Variables
objective = @test_functions;      % Cost Function
 
%% initialize
[row,var] = size(mini);
x = zeros(pop,var);
fnew = zeros(pop,1);
f = zeros(pop,1);
fopt= zeros(pop,1);
xopt=zeros(1,var);

%%  Generation and Initialize the positions 
for i=1:var
    x(:,i) = mini(i)+(maxi(i)-mini(i))*rand(pop,1);
end

for i=1:pop
    f(i) = objective(x(i,:),F_index);
end

%%  Main Loop
gen=1;
while(gen <= maxGen)

    [row,col]=size(x);
    [t,tindex]=min(f);
    Best=x(tindex,:);
    [w,windex]=max(f);
    worst=x(windex,:);
    xnew=zeros(row,col);
    
    for i=1:row
        for j=1:col
            xnew(i,j)=(x(i,j))+rand*(Best(j)-abs(x(i,j))) - (worst(j)-abs(x(i,j)));  % 
        end
    end 
    
    for i=1:row
        xnew(i,:) = max(min(xnew(i,:),maxi),mini);   
        fnew(i,:) = objective(xnew(i,:),F_index);
    end
    
    for i=1:pop
        if(fnew(i)<f(i))
            x(i,:) = xnew(i,:);
            f(i) = fnew(i);
        end
    end

    fnew = []; xnew = [];
    [fopt(gen),ind] = min(f);
    xopt(gen,:)= x(ind,:);
    gen = gen+1;  
    disp(['Iteration No. = ',num2str(gen), ',   Best Cost = ',num2str(min(f))])
    
end

%% 

[val,ind] = min(fopt);
Fes = pop*ind;
disp(['Optimum value = ',num2str(val,10)])



 figure(1)
 plot(fopt,'linewid',2)
 xlabel('Itteration')
 ylabel('Best Cost');
 legend('JAYA')
 disp(' ' )

