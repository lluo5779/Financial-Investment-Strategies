function  x_optimal = MVO_card(mu, Q, targetRet, card, tickers)
    
    % Use this function to construct your MVO portfolio subject to the
    % target return, with short-selling allowed. 
    %
    % You may use Gurobi or any other optimizer that is able to solve mixed
    % integer programs (MIPs). Just be sure to comment on your code to 
    % (briefly) explain your procedure.

    % Find the total number of assets
    n = size(Q,1); 

    % *************** WRITE YOUR CODE HERE ***************
    %----------------------------------------------------------------------
    
    clear model

    names = tickers;

    model.Q = sparse([Q zeros(n,n); zeros(n,2*n)]);
    model.A = sparse([-1*mu', zeros(1,n); 
        ones(1,n), zeros(1,n); 
        zeros(1,n), ones(1,n)]);


    model.obj = [zeros(1,2*n)];
    model.rhs = [-1*targetRet; 1; card];
    sense = [repmat('<', 1, 3)];
    sense(1,[2,3]) = '=';
    model.sense = sense;

    model.vtype = [repmat('C', 1, n) repmat('B', 1, n);];

    gurobi_write(model, 'qp.lp'); % mip.lp

    results = gurobi(model)
% 
%     for v=1:length(names)
%         fprintf('%s %e\n', names{v}, results.x(v));
%     end
% 
%     fprintf('Obj: %e\n', results.objval);



    results  = gurobi(model);
% 
%     for v=1:length(names)
%         fprintf('%s %e\n', names{v}, results.x(v));
%     end

    fprintf('Obj: %e\n', results.objval);


    
    
    
    x_optimal = results.x(1:(size(results.x,1)/2));
    
    %----------------------------------------------------------------------
    
end