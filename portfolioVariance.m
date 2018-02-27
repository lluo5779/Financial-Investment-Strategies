function risk = portfolioVariance(weights, Q)
    n_asset = size(Q,1);
    risk = 0;
    for i = 1:n_asset
        for j = 1:n_asset
            risk = risk + weights(i)*weights(j)*Q(i,j);
        end
    end
    

end