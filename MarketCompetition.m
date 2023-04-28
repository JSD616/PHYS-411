% MarketCompetition.m
% 4-19-23

N_0 = 200; % Initial number of competitors
N_max = 200; % Maximum competitors

active = zeros(1,N_max); % Represents which competitors are active or not

T_max = 30; % Number of iterations
N = zeros(1,T_max); % Number of competitors at each time step
N(1,1) = N_0;

s_a = 5;
mu_a = 40;
s_b = 2;
mu_b = 6;
s_c = 1;
mu_c = -2;
s_d = 5;
mu_d = 10;
startPriceMean = (mu_a - mu_c)./(mu_b + mu_d); % Non-coupled equilibrium price at mean parameter values
%startPriceMean = 500;
startPriceSD = 1;

g = 0.5; % coupling constant (price comparison strength)

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
for n = 1:N_0
    quantity(1,n) = C(1,n) + D(1,n).*startPrice(1,n);
    price(1,n) = (A(1,n)+g.*startAvgPrice-quantity(1,n))./(B(1,n)+g);
    marketValue(1,1) = marketValue(1,1) + price(1,n).*quantity(1,n);
end

for n = 1:N_0
    M(1,n) = price(1,n).*quantity(1,n)./marketValue(1,1);
    avgPrice(1,1) = avgPrice(1,1) + price(1,n).*M(1,n);
end

% Add new competitors randomly
    
% All other steps
for t = 2:T_max
    for n = 1:N_0
        % possible rewards for performance?
        A(t,n) = A(t-1,n);
        B(t,n) = B(t-1,n);
        C(t,n) = C(t-1,n);
        D(t,n) = D(t-1,n);
        quantity(t,n) = C(t,n) + D(t,n).*price(t-1,n);
        price(t,n) = (A(t,n)+g.*avgPrice(t-1,1)-quantity(t))./(B(t,n)+g);
        marketValue(t,1) = marketValue(t,1) + price(t,n).*quantity(t,n);
    end
    
    for n = 1:N_0
        M(t,n) = price(t,n).*quantity(t,n)./marketValue(t,1);
        avgPrice(t,1) = avgPrice(t,1) + price(t,n).*M(t,n);
    end
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

