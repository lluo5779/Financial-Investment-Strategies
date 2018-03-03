function  x_optimal = BL(mu_mkt, Q, tau, P, q, omega, mktNoShares, currentPrices)
    
    % Use this function to construct your MVO portfolio subject to the
    % target return, with short-selling allowed. 
    %
    % You may use quadprog, Gurobi or any other optimizer you are familiar 
    % with. Just be sure to comment on your code to (briefly) explain your 
    % procedure.

    % Find the total number of assets
    n = size(Q,1); 

    % Market portfolio
    x_mkt = mktNoShares .* currentPrices ./ sum( mktNoShares .* currentPrices );
    
    % Variance of market excess returns
    sigma_mkt = x_mkt' * Q * x_mkt;

    % Risk aversion Coefficient
    lambda = mu_mkt/sigma_mkt; % HOW DO I GET THE MARKET EXCESS RETURN???
    
    % Implied market return
    Pi = lambda * Q * x_mkt;
    
    % mu incorporating our views
    u_bar = inv(inv(tau*Q) + P'*inv(omega)*P) * ( inv(tau*Q)*Pi + P'*inv(omega)*q);

    % initializing variables for quadprog           
    f = -1*u_bar;    
    A = [];
    b = [];
    Aeq = [ones(1,n)]; % Budget Constraint
    beq = [1];
    
    [x_optimal, val] = quadprog(0.5*lambda*Q, f, A, b, Aeq, beq, [], [])


    
    
    
    % x_optimal = 
    
    %----------------------------------------------------------------------
    
end