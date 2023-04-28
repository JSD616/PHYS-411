% VarNumMarketComp.m
% 4-23-23

N_0 = 10; % Initial number of competitors
N_max = 50; % Maximum competitors

active = zeros(1,N_max); % Represents which competitors are active or not

T_max = 1000; % Number of iterations
N = zeros(T_max,1); % Number of competitors over time
N(1,1) = N_0;

s_a = 5;
mu_a = 80;
s_b = 2;
mu_b = 16;
s_c = 1;
mu_c = -2;
s_d = 1;
mu_d = 10;
startPriceMean = (mu_a - mu_c)./(mu_b + mu_d); % Non-coupled equilibrium price at mean parameter values
%startPriceMean = 500;
startPriceSD = 1;

g = 2; % price comparison strength parameter
k = 0.09; % growth parameter
p = 1./2.0; % normalization parameter

A = zeros(T_max,N_max);
B = zeros(T_max,N_max);
C = zeros(T_max,N_max);
D = zeros(T_max,N_max);

% For Initializing
startPrice = zeros(1,N_max);
startQuant = zeros(1,N_max);
startMarkVal = 0;
startM = zeros(1,N_max);
startAvgPrice = 0;

price = zeros(T_max,N_max);
quantity  = zeros(T_max,N_max);
marketValue = zeros(T_max,1);
M = zeros(T_max,N_max); % Market shares over time
avgPrice = zeros(T_max,1);

M_eff = zeros(1,N_max); % For normalizing consumer demand

% Initializing step (t = 0)
for n = 1:N_0
    active(1,n) = 1;
    A(1,n) = mu_a + s_a.*randn;
    B(1,n) = mu_b + s_b.*randn;
    C(1,n) = mu_c + s_c.*randn;
    D(1,n) = mu_d + s_d.*randn;
    startPrice(1,n) = startPriceMean + startPriceSD.*randn;
    startQuant(1,n) = C(1,n) + D(1,n).*startPrice(1,n);
    startMarkVal = startMarkVal + startPrice(1,n).*startQuant(1,n);
end

for n = 1:N_0
    startM(1,n) = startPrice(1,n).*startQuant(1,n)./startMarkVal;
    startAvgPrice = startAvgPrice + startPrice(1,n).*startM(1,n);
end

% First step (t = 1)

% Effective market shares to normalize consumer demand
for n = 1:N_0
    M_eff(1,n) = max([startM(1,n) 1./(2.*N_max)]);
    %M_eff(1,n) = 1;
end
M_eff = M_eff.^p;
for n = 1:N_0
    quantity(1,n) = C(1,n) + D(1,n).*startPrice(1,n);
    price(1,n) = (A(1,n)+g.*startAvgPrice-quantity(1,n)./M_eff(1,n))./(B(1,n)+g);
    
    % Competitors may exit
    if (active(1,n)==1) && ((quantity(1,n)<0) || (price(1,n)<0))
        quantity(1,n) = 0;
        price(1,n) = 0;
        active(1,n) = 0;
        N(1,1) = N(1,1) - 1;
        %N(1,1)
        %n
    end
    marketValue(1,1) = marketValue(1,1) + price(1,n).*quantity(1,n);
end

for n = 1:N_0
    M(1,n) = price(1,n).*quantity(1,n)./marketValue(1,1);
    avgPrice(1,1) = avgPrice(1,1) + price(1,n).*M(1,n);
end

% Add new competitors randomly

lambda = log(marketValue(1,1).*N(1,1) + 1).*k;
%lambda
delta_N = poissrnd(lambda);
delta_N = min([delta_N N_max - N(1,1)]);
%delta_N
N(1,1) = N(1,1) + delta_N;

n = 1;

while (delta_N > 0)
    %delta_N
    if (active(1,n)==0)
        %n
        A(1,n) = mu_a + s_a.*randn;
        B(1,n) = mu_b + s_b.*randn;
        C(1,n) = mu_c + s_c.*randn;
        D(1,n) = mu_d + s_d.*randn;
        price(1,n) = avgPrice(1,1);
        M(1,n) = 0;
        active(1,n) = 1;
        delta_N = delta_N - 1;
        %delta_N
    end
    n = n + 1;    
end
% All other steps
for t = 2:T_max
    N(t,1) = N(t-1,1);
    
    % Effective market shares to normalize consumer demand
    for n = 1:N_max
        M_eff(1,n) = max([M(t-1,n) 1./(2.*N_max)]);
        %M_eff(1,n) = 1;
    end
    
    M_eff = M_eff.^p;
    
    for n = 1:N_max
        % possible rewards for performance?
        A(t,n) = A(t-1,n);
        B(t,n) = B(t-1,n);
        C(t,n) = C(t-1,n);
        D(t,n) = D(t-1,n);
        if (active(1,n)==1)
            quantity(t,n) = C(t,n) + D(t,n).*price(t-1,n);
            price(t,n) = (A(t,n)+g.*avgPrice(t-1,1)-quantity(t)./M_eff(1,n))./(B(t,n)+g);
        end
        % Competitors may exit
        if (active(1,n)==1) && ((quantity(t,n)<0) || (price(t,n)<0))
            %n
            quantity(t,n) = 0;
            price(t,n) = 0;
            active(1,n) = 0;
            N(t,1) = N(t,1) - 1;
        end
        marketValue(t,1) = marketValue(t,1) + price(t,n).*quantity(t,n);
    end
    
    for n = 1:N_max
        M(t,n) = price(t,n).*quantity(t,n)./marketValue(t,1);
        avgPrice(t,1) = avgPrice(t,1) + price(t,n).*M(t,n);
    end
    
    % Add new competitors randomly
    lambda = log(marketValue(t,1).*N(t,1) + 1).*k;
    delta_N = poissrnd(lambda);
    delta_N = min([delta_N N_max-N(t,1)]);
    %delta_N
    N(t,1) = N(t,1) + delta_N;
    %N(t,1)
    
    n = 1;
    
    while (delta_N > 0)
        if (active(1,n)==0)
            %n
            A(t,n) = mu_a + s_a.*randn;
            B(t,n) = mu_b + s_b.*randn;
            C(t,n) = mu_c + s_c.*randn;
            D(t,n) = mu_d + s_d.*randn;
            price(t,n) = avgPrice(t,1);
            M(t,n) = 0;
            active(1,n) = 1;
            delta_N = delta_N - 1;
        end
        n = n + 1;    
    end
    
    for n = 1:N_max
        if (M(t,n) < 0)
            n
        end
    end
    
    %disp ('--------');
end

figure(1)
plot(price)
title('Price')

figure(2)
plot(avgPrice)
title('Average Price')

figure(3)
plot(quantity)
title('Quantity Produced')

figure(4)
plot(marketValue)
title('Market Value')

figure(5)
plot(M)
title('Market Share')

figure(6)
plot(N)
title('Number of Competitors')
