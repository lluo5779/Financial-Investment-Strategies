    %% MIE377 (Winter 2018) - Project 1
% The purpose of this program is to test the out-of-sample performance of 
% 3 different portfolio optimization models. We will test the following
% models:
% 1. MVO
% 2. MVO with cardinality constraint
% 3. Black-Litterman model
%
% Use this template to write your program
%
% Student Name: Yiqing (Louis) Luo
% Student ID:   1002449059

clc
clear all
format long

% Program Start
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Read input files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load the stock weekly prices and factors weekly returns
adjClose = readtable('Project1_Data_adjClose.csv');%, 'ReadRowNames', true);
adjClose.Properties.RowNames = cellstr(datetime(adjClose.Date));
%adjClose.Properties.RowNames = datetime(cellstr(table2array(adjClose(:,1)))); %cellstr(datetime(adjClose.Properties.RowNames));
size_adjClose = size(adjClose);
adjClose = adjClose(:,2:size_adjClose(2));

factorRet = readtable('Project1_Data_FF_factors.csv'); %, 'ReadRowNames', true);
factorRet.Properties.RowNames = cellstr(datetime(factorRet.Date));
size_factorRet = size(factorRet);
factorRet = factorRet(:,2:size_factorRet(2));

riskFree = factorRet(:,4);
factorRet = factorRet(:,1:3);

% Identify the tickers and the dates 
tickers = adjClose.Properties.VariableNames';
dates   = datetime(factorRet.Properties.RowNames);

% Calculate the stocks' weekly EXCESS returns
prices  = table2array(adjClose);
returns_raw = ( prices(2:end,:) - prices(1:end-1,:) ) ./ prices(1:end-1,:);
returns = returns_raw - ( diag( table2array(riskFree) ) * ones( size(returns_raw) ) );
returns = array2table(returns);
returns.Properties.VariableNames = tickers;
returns.Properties.RowNames = cellstr(datetime(factorRet.Properties.RowNames));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Define your initial parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of shares outstanding per company
mktNoShares = [36.40 9.86 15.05 10.91 45.44 15.43 43.45 33.26 51.57 45.25 ...
                69.38 8.58 64.56 30.80 36.92 7.70 3.32 54.68 38.99 5.78]' * 1e8;

% Initial budget to invest
initialVal = 100;

% Start of in-sample calibration period 
calStart = datetime('2012-01-01');
calEnd   = calStart + calmonths(12) - days(1);

% Start of out-of-sample test period 
testStart = datetime('2013-01-01');
testEnd   = testStart + calmonths(6) - days(1);
       
% Number of investment periods (each investment period is 6 months long)
NoPeriods = 6;

% Investment strategies
% Note: You must populate the functios MVO.m, MVO_card.m and BL.m with your
% own code to construct the optimal portfolios. 
funNames  = {'MVO' 'MVO (Card=12)' 'B-L Model'};
NoMethods = length(funNames);

funList = {'MVO' 'MVO_card' 'BL'};
funList = cellfun(@str2func, funList, 'UniformOutput', false);

% Maximum number of assets imposed by cardinality constraint
card = 12;

% *************** WRITE YOUR CODE HERE ***************
%--------------------------------------------------------------------------

% Basic visualization of return over initial training period
NoShares{1} = mktNoShares;
NoShares{2} = mktNoShares;
NoShares{3} = mktNoShares;

[NoTotalDates, NoAssets] = size(adjClose);
% Insert your B-L parameters (have at least 6 different views):
% 
% tau = 0.005;
% P = ones(1,NoAssets)*1/NoAssets;
% q = [geomean(1+geomean(1+periodReturns(NoTotalDates/2:end,NoAssets))-1)-1]
for i = 1:NoMethods
    avgPortfolioValue{i} = zeros(NoPeriods,1);
    expectedPortfReturn{i} = zeros(NoPeriods,1);
end

%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Construct and rebalance your portfolios
%
% Here you will estimate your input parameters (exp. returns, cov. matrix,
% etc) from the Fama-French factor models. You will have to re-estimate 
% your parameters at the start of each rebalance period, and then 
% re-optimize and rebalance your portfolios accordingly. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initiate counter for the number of observations per investment period
toDay = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CODE CHANGED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NoMethods = 3;
currentVal = ones(NoPeriods, NoMethods);

NoAbsConstraints = 3;
NoRelaConstraints = 5;

% indicesSimilarCompanies{1} = [1,2]
% indicesSimilarCompanies{2} = [4,5,6]
% indicesSimilarCompanies{3} = [5,6]
% indicesSimilarCompanies{4} = [8,9,10]
% indicesSimilarCompanies{5} = [8,9]
% indicesSimilarCompanies{6} = [9,10]
% indicesSimilarCompanies{7} = [8,10]
% indicesSimilarCompanies{8} = [11,12]
% indicesSimilarCompanies{9} = [13,14]
% indicesSimilarCompanies{10} = [15,16,17,20]
% indicesSimilarCompanies{11} = [15,16,17]
% indicesSimilarCompanies{12} = [15,16,20]
% indicesSimilarCompanies{13} = [15,17, 20]
% indicesSimilarCompanies{13} = [16,17, 20]
% indicesSimilarCompanies{14} = [15,16]
% indicesSimilarCompanies{15} = [15,17]
% indicesSimilarCompanies{16} = [15, 20]
% indicesSimilarCompanies{17} = [15,20]
% indicesSimilarCompanies{18} = [16,17]
% indicesSimilarCompanies{19} = [16,20]
% indicesSimilarCompanies{20} = [17,20]
% indicesSimilarCompanies{21} = [19,20]

indicesSimilarCompanies = nchoosek(linspace(1,NoAssets, NoAssets),2);

for t = 1 : NoPeriods
    
    % Subset the returns and factor returns corresponding to the current
    % calibration period.
    periodReturns = table2array( returns( calStart <= dates & dates <= calEnd, :) );
    periodFactRet = table2array( factorRet( calStart <= dates & dates <= calEnd, :) );
    currentPrices = table2array( adjClose( ( calEnd - days(7) ) <= dates ... 
                                                    & dates <= calEnd, :) )';
    currentReturns = table2array( returns( ( calEnd - days(7) ) <= dates ... 
                                                    & dates <= calEnd, :));
    startingReturns = table2array( returns( calStart <= dates ... 
                                                    & dates <= (calStart +   days(7)), :));
    avgPeriodRiskFree = table2array(riskFree(( calEnd - days(7) ) <= dates ... 
                                                    & dates <= calEnd, :));
    currentRiskFree = table2array( riskFree( ( calEnd - days(7) ) <= dates ... 
                                                    & dates <= calEnd, :));
    startingPrices = table2array( adjClose( calStart <= dates ... 
                                                    & dates <= (calStart +   days(7)), :))';                                           
    % Subset the prices corresponding to the current out-of-sample test 
    % period.
    periodPrices = table2array( adjClose( testStart <= dates & dates <= testEnd,:) );
    
    % Set the initial value of the portfolio or update the portfolio value
    if t == 1
        
        currentVal(t,:) = initialVal;
        
    else
        for i = 1 : NoMethods
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CODE CHANGED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %currentVal(t,i) = currentPrices .* NoShares{i};
            currentVal(t,i) = currentPrices' * NoShares{i}; %<=CHANGED GOOD CHANGE(Y)
            
        end
    end

    
    
    % Update counter for the number of observations per investment period
    fromDay = toDay + 1;
    toDay   = toDay + size(periodPrices,1);
    
    % *************** WRITE YOUR CODE HERE ***************
    %----------------------------------------------------------------------
    
    % Calculate your initial exp. returns and covariance matrix using the
    % Fama-French factor model. You may find these values as you prefer,
    % the code given below is just an example. 
    
    % n = 20 assets
    % m = 3
    % for each asset, have a F-F model representation
    
    % A =           % n x 1 vector of alphas
    % V =           % m x n matrix of betas
    % 
    % f_bar =       % m x 1 vector of factor expected returns
    % F =           % m x m factor covariance matrix
    % 
    % epsilon =     % Regression residuals
    % D =           % Diagonal n x n matrix of residual variance

    % mu =          % n x 1 vector of asset exp. returns
    % Q  =          % n x n asset covariance matrix
    
    
    X = [ones(size(periodFactRet,1),1) periodFactRet];
    B = (X'*X)\(X'*periodReturns); %inv(X'*X)*X'*periodReturns
    A = B(1,:)'; %[b0; b0; b0;...]
    V = B(2:size(B,1),:)'; %[b1 b2 b3; b1 b2 b3; ...]
    
    f_bar = [1 (geomean(X(:,2:end)+1) - 1)]'; % [f_bar1; f_bar2; f_bar3]
    F = cov(periodFactRet);
    predicted_ret = X * B;
    
    epsilon = predicted_ret - periodReturns;
    om = cov(epsilon);
    D = diag(diag(om));
    
    mu = B'*f_bar; %[u1;u2;u3; ...]
    Q = V*F*V' + om;
  
%     var_mkt = var(table2array(factorRet(:,1)));
%     lambda = mu_mkt/var_mkt;




    %------------ Black Litterman Parameters -------------%
      
    tau = 0.001;
    P = zeros(NoAbsConstraints + NoRelaConstraints, NoAssets);
    q = zeros(NoAbsConstraints + NoRelaConstraints, 1);
   
    mu_mkt = geomean(table2array(factorRet(:,1))+1)-1;
    
    [mu_sorted1, I1] = sort(mu, 'descend');
    
    for i = 1:NoAbsConstraints
        P(i, I1(i,1)) = 1;
        q(i,1) = 2*(mu_sorted1(i,1));
    end
    
    for i = 1:size(indicesSimilarCompanies,1)
        mu_1 = mu(indicesSimilarCompanies(i,1),1);
        mu_2 = mu(indicesSimilarCompanies(i,2),1);
        abs_diff_mu(i,:) = abs(mu_1 - mu_2);
        diff_mu(i,:) = mu_1 - mu_2;
    end
    
    [mu_sorted2, I2] = sort(diff_mu, 'descend');
    
    for i = 1:NoRelaConstraints
        P(NoAbsConstraints+i, indicesSimilarCompanies(I2(i,:),1)) = 1;
        P(NoAbsConstraints+i, indicesSimilarCompanies(I2(i,:),2)) = -1;
        q(NoAbsConstraints+i,1) = diff_mu(I2(i,1),1)*2;
    end
    
    Omega = diag(diag(tau*P*Q*P'));
        
    %----------------------------------------------------------------------
    
    % Define the target return for the 2 MVO portfolios
    targetRet = mean(mu);
    
    % Optimize your portfolios to get the weights 'x'
    % Note: You need to write the code for the 3 functions ('MVO.m', 
    % 'MVO_card.m' and 'BL.m') to find your optimal portfolios
    
    x{1}(:,t) = funList{1}(mu, Q, targetRet); 
    x{2}(:,t) = funList{2}(mu, Q, targetRet, card, tickers); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CODE CHANGED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x{3}(:,t) = funList{3}(mu_mkt, Q, tau, P, q, Omega, mktNoShares, currentPrices); 
    
    
    % Calculate the optimal number of shares of each stock you should hold
    for i = 1:NoMethods
        
        % Number of shares your portfolio holds per stock
        NoShares{i} = x{i}(:,t) .* currentVal(t,i) ./ currentPrices;
        
        % Weekly portfolio value during the out-of-sample window
        portfValue(fromDay:toDay,i) = periodPrices * NoShares{i};
        
        % *************** WRITE YOUR CODE HERE ***************
        %------------------------------------------------------------------
        
        % Calculate your transaction costs for the current rebalance
        % period. The first period does not have any cost since you are
        % constructing the portfolios for the 1st time. 
        
        if t ~= 1
           
           abs_diff = abs(NoSharesOld{i} - NoShares{i});
           tCost(t-1, i) = abs_diff' * currentPrices * 0.005; % 0.5% of traded volume
            
        end
        
        NoSharesOld{i} = NoShares{i};
        
        % The average portfolio value is determined by finding the
        % arithmetic mean of all the values of portfolio during the period
        avgPortfolioValue{i}(t,1) = mean(portfValue(fromDay:toDay,i));
        % Calculation of Sharpe Ratio - the difference between the expected 
        % values of the market return and risk free rate in the week prior
        % divided between the portfolio variance over that period of time.
        expectedCurrentMarketReturn = currentReturns*x{i}(:,t); %geomean(1+currentReturns)-1;
        mu_portfolio{i}(t,:) = mu'*x{i}(:,t);
        portfolio_var{i}(t,:) = x{i}(:,t)'*Q*x{i}(:,t);
        %portfolio_var_sortino{i}(t,:) = x{i}(:,t)'*cov(min(predicted_ret,0))*x{i}(:,t);
        sharpe_ratio{i}(t,:) = (mu_portfolio{i}(t,:))/sqrt(portfolio_var{i}(t,:));
        %sortino_ratio{i}(t,:) = (mu_portfolio{i}(t,:))/sqrt(portfolio_var_sortino{i}(t,:));
        
        
    end

    % Update your calibration and out-of-sample test periods
    calStart = calStart + calmonths(6);
    calEnd   = calStart + calmonths(12) - days(1);
    
    testStart = testStart + calmonths(6);
    testEnd   = testStart + calmonths(6) - days(1);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% *************** WRITE YOUR CODE HERE ***************
%--------------------------------------------------------------------------

% Calculate the portfolio average return, variance (or standard deviation),
% or any other performance and/or risk metric you wish to include in your
% report.

%--------------------------------------------------------------------------
% for i = 1:NoMethods
%     avgPortfolioValue{i}
%     sharpe_ratio{i}
%     
% end   
        

save('portfolioParameters.mat', 'avgPortfolioValue', 'sharpe_ratio', 'expectedCurrentMarketReturn', 'mu_portfolio','portfolio_var','-v7.3')
% save('sharpe_ratio.mat', sharpe_ratio)
% save('expectedCurrentMarketReturn.mat',expectedCurrentMarketReturn)
% save('expectedmu.mat',expectedmu)

%--------------------------------------------------------------------------
% 4.1 Plot the portfolio values 
%--------------------------------------------------------------------------
testStart = datetime('2013-01-01');
plotDates = dates(testStart <= dates);

fig1 = figure(1);
plot(plotDates, portfValue(:,1))
hold on
plot(plotDates, portfValue(:,2))
hold on
plot(plotDates, portfValue(:,3))
legend(funNames, 'Location', 'eastoutside','FontSize',12);
datetick('x','dd-mmm-yyyy','keepticks','keeplimits');
set(gca,'XTickLabelRotation',30);
title('Portfolio Value', 'FontSize', 14)
ylabel('Value','interpreter','latex','FontSize',12);

% Define the plot size in inches
set(fig1,'Units','Inches', 'Position', [0 0 8, 5]);
pos1 = get(fig1,'Position');
set(fig1,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos1(3), pos1(4)]);

% If you want to save the figure as .pdf for use in LaTeX
% print(fig1,'fileName','-dpdf','-r0');

% If you want to save the figure as .png for use in MS Word
print(fig1,'protfValue','-dpng','-r0');

%--------------------------------------------------------------------------
% 4.2 Plot the portfolio weights 
%--------------------------------------------------------------------------

% MVO Plot
fig2 = figure(2);
area(x{1}')
legend(tickers, 'Location', 'eastoutside','FontSize',12);
title('Portfolio Weights', 'FontSize', 14)
ylabel('Weights','interpreter','latex','FontSize',12);
xlabel('Rebalance Period','interpreter','latex','FontSize',12);

% Define the plot size in inches
set(fig2,'Units','Inches', 'Position', [0 0 8, 5]);
pos1 = get(fig2,'Position');
set(fig2,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos1(3), pos1(4)]);

% If you want to save the figure as .pdf for use in LaTeX
% print(fig2,'fileName2','-dpdf','-r0');

% If you want to save the figure as .png for use in MS Word
print(fig2,'weights_mvo','-dpng','-r0');


% MVO with Cardinality Constraints Plot
fig3 = figure(3);
area(x{2}')
legend(tickers, 'Location', 'eastoutside','FontSize',12);
title('Portfolio Weights', 'FontSize', 14)
ylabel('Weights','interpreter','latex','FontSize',12);
xlabel('Rebalance Period','interpreter','latex','FontSize',12);

% Define the plot size in inches
set(fig3,'Units','Inches', 'Position', [0 0 8, 5]);
pos1 = get(fig3,'Position');
set(fig3,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos1(3), pos1(4)]);

% If you want to save the figure as .pdf for use in LaTeX
% print(fig3,'fileName3','-dpdf','-r0');

% If you want to save the figure as .png for use in MS Word
print(fig3,'weights_card','-dpng','-r0');


% B-L Model Plot
fig4 = figure(4);
area(x{3}')
legend(tickers, 'Location', 'eastoutside','FontSize',12);
title('Portfolio Weights', 'FontSize', 14)
ylabel('Weights','interpreter','latex','FontSize',12);
xlabel('Rebalance Period','interpreter','latex','FontSize',12);

% Define the plot size in inches
set(fig4,'Units','Inches', 'Position', [0 -5 8, 5]);
pos1 = get(fig4,'Position');
set(fig4,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos1(3), pos1(4)]);

% If you want to save the figure as .pdf for use in LaTeX
%print(fig4,'fileName3','-dpdf','-r0');

% If you want to save the figure as .png for use in MS Word
print(fig4,'weights_BL','-dpng','-r0');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MVO Plot
fig5 = figure(5);
plot(x{1}')
legend(tickers, 'Location', 'eastoutside','FontSize',12);
title('Portfolio Weights', 'FontSize', 14)
ylabel('Weights','interpreter','latex','FontSize',12);
xlabel('Rebalance Period','interpreter','latex','FontSize',12);

% Define the plot size in inches
set(fig5,'Units','Inches', 'Position', [0 0 8, 5]);
pos1 = get(fig5,'Position');
set(fig5,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos1(3), pos1(4)]);

% If you want to save the figure as .pdf for use in LaTeX
% print(fig2,'fileName2','-dpdf','-r0');

% If you want to save the figure as .png for use in MS Word
print(fig5,'weights_mvo2','-dpng','-r0');


% MVO with Cardinality Constraints Plot
fig6 = figure(6);
plot(x{2}')
legend(tickers, 'Location', 'eastoutside','FontSize',12);
title('Portfolio Weights', 'FontSize', 14)
ylabel('Weights','interpreter','latex','FontSize',12);
xlabel('Rebalance Period','interpreter','latex','FontSize',12);

% Define the plot size in inches
set(fig6,'Units','Inches', 'Position', [0 0 8, 5]);
pos1 = get(fig6,'Position');
set(fig6,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos1(3), pos1(4)]);

% If you want to save the figure as .pdf for use in LaTeX
% print(fig3,'fileName3','-dpdf','-r0');

% If you want to save the figure as .png for use in MS Word
print(fig6,'weights_card2','-dpng','-r0');


% B-L Model Plot
fig7 = figure(7);
plot(x{3}')
legend(tickers, 'Location', 'eastoutside','FontSize',12);
title('Portfolio Weights', 'FontSize', 14)
ylabel('Weights','interpreter','latex','FontSize',12);
xlabel('Rebalance Period','interpreter','latex','FontSize',12);

% Define the plot size in inches
set(fig7,'Units','Inches', 'Position', [0 -5 8, 5]);
pos1 = get(fig7,'Position');
set(fig7,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos1(3), pos1(4)]);

% If you want to save the figure as .pdf for use in LaTeX
%print(fig4,'fileName3','-dpdf','-r0');

% If you want to save the figure as .png for use in MS Word
print(fig7,'weights_BL2','-dpng','-r0');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program End