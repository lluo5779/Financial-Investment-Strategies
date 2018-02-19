function  x_optimal = BL(Q, tau, P, q, Omega, mktNoShares, currentPrices)
    
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
    
    % *************** WRITE YOUR CODE HERE ***************
    %----------------------------------------------------------------------
    
    
    
    
    
    
    
    % x_optimal = 
    
    %----------------------------------------------------------------------
    
end