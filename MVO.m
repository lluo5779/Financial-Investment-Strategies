function  x_optimal = MVO(mu, Q, targetRet)
    
    % Use this function to construct your MVO portfolio subject to the
    % target return, with short-selling allowed. 
    %
    % You may use quadprog, Gurobi, or any other optimizer you are familiar
    % with. Just be sure to comment on your code to (briefly) explain your
    % procedure.

    % Find the total number of assets
    n = size(Q,1); 

    % *************** WRITE YOUR CODE HERE ***************
    %---------------------------------------------------------------------
    
    A = [];
    Aeq = [-1*mu';ones(1,n)];
    
    b = [];
    beq = [-1*targetRet;1];
    
    
    [x_optimal, val] = quadprog(Q, [], A, b, Aeq, beq, [], []);
    
    
    
    %----------------------------------------------------------------------
    
end