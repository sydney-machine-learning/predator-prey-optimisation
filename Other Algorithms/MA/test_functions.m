% original
function fit=test_functions(L,F_index)
%Insert your own objective function with a new F_index.

if F_index==1
    dim = length(L);
    term1 = sum(L.^2);
    term2 = (0.5 * sum((1:dim) .* L))^2;
    term3 = (0.5 * sum((1:dim) .* L))^4;
    fit = term1 + term2 + term3;
end

if F_index==2 
    fit = 0;
    dim = length(L);
    for i = 1:(dim-1)
        term = 100 * (L(i+1) - L(i)^2)^2 + (L(i) - 1)^2;
        fit = fit + term;
    end
end

if F_index==3
    dim = length(L);
    fit = 0;
    for i = 1:(dim-1)
        term = 0.5 + ((sin(sqrt(L(i)^2 + L(i+1)^2))^2 - 0.5) / (1 + 0.001 * (L(i)^2 + L(i+1)^2))^2);
        fit = fit + term;
    end
end

if F_index==4
    dim = length(L);
    term1 = 10 * dim;
    term2 = sum(L.^2 - 10 * cos(2 * pi * L));
    fit = term1 + term2;
end

if F_index==5
    dim = length(L);
    L=1+(L-1)/4;
    term1 = (sin(pi * L(1)))^2;
    term2 = sum(((L(1:dim-1)-1).^2).*(1 + 10 * (sin(pi * L(1:dim-1) - 1)).^2));
    term3 = (L(dim - 1) - 1)^2;
    term4 = 1 + sin(2 * pi * L(dim))^2;
    
    fit = term1 + term2 + term3 * term4;
end

if F_index==6 
    D = length(L);
    fit = L(1)^2 + 1e6 * sum(L(2:D).^2);
end

if F_index==7
    D = length(L);
    fit=(sum(L.^2)^2-sum(L)^2)^0.5 + (0.5*sum(L.^2)+sum(L))/D + 0.5; 
end

if F_index==8
    D = length(L);
    fit = 0;
    for i = 1: D
        fit = fit +(L(i)^2)*(1e6)^((i-1)/(D-1));
    end
end

if F_index==9
    D= length(L);
    fit=(sum(L.^2)^2-D)^0.25 + (0.5*sum(L.^2)+sum(L))/D + 0.5; 
end

if F_index==10
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

if F_index==11
    % Ackley function
    D = length(L);
    sum1 = sum(L.^2);
    sum2 = sum(cos(2 * pi * L));
    term1 = -20 * exp(-0.2 * sqrt(sum1 / D));
    term2 = -exp(sum2 / D);
    fit = term1 + term2 + 20 + exp(1);
end

if F_index==12
    % Discus function
    sum_term = sum(L(2:end).^2);
    fit = 1e6 * L(1)^2 + sum_term;
end

if F_index==13
    % Griewank function
    D = length(L);
    sum_term = sum(L.^2) / 4000;
    prod_term = prod(cos(L ./ sqrt(1:D)));
    fit = sum_term - prod_term + 1;
end

if F_index==14
    % Schaffer F14 function
    D = length(L);
    sum_term = 0;
    for i = 1:(D-1)
        si = sqrt(L(i)^2 + L(i+1)^2);
        sum_term = sum_term + (sqrt(si) *( sin(50.0 * si^0.2) + 1))^2;
    end
    fit=sum_term/(D-1);
end
end