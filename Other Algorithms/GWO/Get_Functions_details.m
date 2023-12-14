%  Sine Cosine Algorithm (SCA)  
%
%  Source codes demo version 1.0                                                                      
%                                                                                                     
%  Developed in MATLAB R2011b(7.13)                                                                   
%                                                                                                     
%  Author and programmer: Seyedali Mirjalili                                                          
%                                                                                                     
%         e-Mail: ali.mirjalili@gmail.com                                                             
%                 seyedali.mirjalili@griffithuni.edu.au                                               
%                                                                                                     
%       Homepage: http://www.alimirjalili.com                                                         
%                                                                                                     
%  Main paper:                                                                                        
%  S. Mirjalili, SCA: A Sine Cosine Algorithm for solving optimization problems
%  Knowledge-Based Systems, DOI: http://dx.doi.org/10.1016/j.knosys.2015.12.022

% This function containts full information and implementations of the benchmark 
% functions in Table 1, Table 2, and other test functins from the literature 

% lb is the lower bound: lb=[lb_1,lb_2,...,lb_d]
% up is the uppper bound: ub=[ub_1,ub_2,...,ub_d]
% dim is the number of variables (dimension of the problem)

function [lb,ub,dim,fobj] = Get_Functions_details(F)


switch F
    case 'F1'
        fobj = @F1;
        lb=-100;
        ub=100;
        dim=100;
        
    case 'F2'
        fobj = @F2;
        lb=-10;
        ub=10;
        dim=100;
        
    case 'F3'
        fobj = @F3;
        lb=-100;
        ub=100;
        dim=100;
        
    case 'F4'
        fobj = @F4;
        lb=-100;
        ub=100;
        dim=100;
        
    case 'F5'
        fobj = @F5;
        lb=-30;
        ub=30;
        dim=100;
        
    case 'F6'
        fobj = @F6;
        lb=-100;
        ub=100;
        dim=100;
        
    case 'F7'
        fobj = @F7;
        lb=-1.28;
        ub=1.28;
        dim=100;
        
    case 'F8'
        fobj = @F8;
        lb=-500;
        ub=500;
        dim=100;
        
    case 'F9'
        fobj = @F9;
        lb=-32;
        ub=32;
        dim=100;

    case 'F10'
        fobj = @F10;
        lb=-50;
        ub=50;
        dim=100;
        
    case 'F11'
        fobj = @F11;
        lb=-50;
        ub=50;
        dim=100;
        
    case 'F12'
        fobj = @F12;
        lb=-65.536;
        ub=65.536;
        dim=100;
        
    case 'F13'
        fobj = @F13;
        lb=-5;
        ub=5;
        dim=100;
        
    case 'F14'
        fobj = @F14;
        lb=-5;
        ub=5;
        dim=100;
        
end

end

% F1

function fit = F1(L)
dim = length(L);
    term1 = sum(L.^2);
    term2 = (0.5 * sum((1:dim) .* L))^2;
    term3 = (0.5 * sum((1:dim) .* L))^4;
    fit = term1 + term2 + term3;
end

% F2

function fit = F2(L)
fit = 0;
    dim = length(L);
    for i = 1:(dim-1)
        term = 100 * (L(i+1) - L(i)^2)^2 + (L(i) - 1)^2;
        fit = fit + term;
    end
end

% F3

function fit = F3(L)
dim = length(L);
    fit = 0;
    for i = 1:(dim-1)
        term = 0.5 + ((sin(sqrt(L(i)^2 + L(i+1)^2))^2 - 0.5) / (1 + 0.001 * (L(i)^2 + L(i+1)^2))^2);
        fit = fit + term;
    end
end

% F4

function fit = F4(L)
dim = length(L);
    term1 = 10 * dim;
    term2 = sum(L.^2 - 10 * cos(2 * pi * L));
    fit = term1 + term2;
end

% F5

function fit = F5(L)
dim = length(L);
    L=1+(L-1)/4;
    term1 = (sin(pi * L(1)))^2;
    term2 = sum(((L(1:dim-1)-1).^2).*(1 + 10 * (sin(pi * L(1:dim-1) - 1)).^2));
    term3 = (L(dim - 1) - 1)^2;
    term4 = 1 + sin(2 * pi * L(dim))^2;
    
    fit = term1 + term2 + term3 * term4;
end

% F6

function fit = F6(L)
D = length(L);
    fit = L(1)^2 + 1e6 * sum(L(2:D).^2);
end

% F7

function fit = F7(L)
D = length(L);
    fit=(sum(L.^2)^2-sum(L)^2)^0.5 + (0.5*sum(L.^2)+sum(L))/D + 0.5; 
end

% F8

function fit = F8(L)
D = length(L);
    fit = 0;
    for i = 1: D
        fit = fit +(L(i)^2)*(1e6)^((i-1)/(D-1));
    end
end

% F9

function fit = F9(L)
D= length(L);
    fit=(sum(L.^2)^2-D)^0.25 + (0.5*sum(L.^2)+sum(L))/D + 0.5; 
end

% F10

function fit = F10(L)
D = length(L);
    % Shift the input vector
    z = L + 420.9687462275036;
    % Calculate the Modified Schwefel's function
    sum_term = 0;
    for i = 1:D
        zi = z(i);
        if abs(zi) <= 500
            sum_term = sum_term + zi * sin(sqrt(abs(zi)));
        elseif zi > 500
            sum_term = sum_term + (500 - mod(zi, 500)) * sin(sqrt(abs(500 - mod(zi, 500)))) - (zi - 500)^2 / 10000*D;
        else
            sum_term = sum_term + (mod(abs(zi), 500) - 500) * sin(sqrt(abs(mod(zi, 500) - 500))) - (zi + 500)^2 / 10000*D;
        end
    end
    % Final objective function value
    fit = 418.9829*D - sum_term;
end

% F11

function fit = F11(L)
% Ackley function
    D = length(L);
    sum1 = sum(L.^2);
    sum2 = sum(cos(2 * pi * L));
    term1 = -20 * exp(-0.2 * sqrt(sum1 / D));
    term2 = -exp(sum2 / D);
    fit = term1 + term2 + 20 + exp(1);
end

% F12

function fit = F12(L)
% Discus function
    sum_term = sum(L(2:end).^2);
    fit = 1e6 * L(1)^2 + sum_term;
end

% F13

function fit = F13(L)
% Griewank function
    D = length(L);
    sum_term = sum(L.^2) / 4000;
    prod_term = prod(cos(L ./ sqrt(1:D)));
    fit = sum_term - prod_term + 1;
end

% F14

function fit = F14(L)
% Schaffer F14 function
    D = length(L);
    sum_term = 0;
    for i = 1:(D-1)
        si = sqrt(L(i)^2 + L(i+1)^2);
        sum_term = sum_term + (sqrt(si) *( sin(50.0 * si^0.2) + 1))^2;
    end
    fit=sum_term/(D-1);
end
