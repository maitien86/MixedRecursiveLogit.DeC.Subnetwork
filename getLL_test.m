% Compute the loglikelohood value and its gradient.
%%
function [LL, grad] = getLL_test()

    global incidenceFull; 
    global Gradient;
    global Op;
    global Mfull;
    global Ufull;
    global Atts;
    global Obs;     % Observation
    global nbobs;  
    global isLinkSizeInclusive;
    global lastIndexNetworkState;
    % ----------------------------------------------------
    % If Link size is included
    mu = 1; 
    % MU IS NORMALIZED TO ONE
    [lastIndexNetworkState, maxDest] = size(incidenceFull);
    
    Mfull = getM(Op.x, isLinkSizeInclusive);
    MregularNetwork = Mfull(1:lastIndexNetworkState,1:lastIndexNetworkState);
    Ufull = getU(Op.x, isLinkSizeInclusive);
    UregularNetwork = Ufull(1:lastIndexNetworkState,1:lastIndexNetworkState);
    % Set LL value
    LL = 0;
    grad = zeros(1, Op.n);
    % Initialize
    M = MregularNetwork;
    M(:,lastIndexNetworkState+1) = sparse(zeros(lastIndexNetworkState,1));
    M(lastIndexNetworkState+1,:) = sparse(zeros(1, lastIndexNetworkState + 1));
    U = UregularNetwork;
    U(:,lastIndexNetworkState+1) = sparse(zeros(lastIndexNetworkState,1));
    U(lastIndexNetworkState+1,:) = sparse(zeros(1, lastIndexNetworkState + 1));
    for i = 1:Op.n
        AttLc(i) =  Matrix2D(Atts(i).Value(1:lastIndexNetworkState,1:lastIndexNetworkState));
        AttLc(i).Value(:,lastIndexNetworkState+1) = sparse(zeros(lastIndexNetworkState,1));
        AttLc(i).Value(lastIndexNetworkState+1,:) = sparse(zeros(1, lastIndexNetworkState + 1));
    end

    % b matrix:
    N = size(M,1);
    b = sparse(zeros(N,1));
    b(N) = 1;
    B = sparse(zeros(N, maxDest - lastIndexNetworkState));
    B(N,:) = ones(1,maxDest - lastIndexNetworkState);
    for i = 1: maxDest - lastIndexNetworkState
        B(1:lastIndexNetworkState,i) = Mfull(:, i+lastIndexNetworkState);
    end
    A = speye(size(M)) - M;
    Z = A\B;
    % Check feasible
    minele = min(Z(:));
    expVokBool = 1;
    if minele == 0 || minele < OptimizeConstant.NUM_ERROR
       expVokBool = 0;
    end 
    Zabs = abs(Z); % MAYBE SET TO VERY SMALL POSITIVE VALUE? 
    D = (A * Z - B);
    resNorm = norm(D(:));
    if resNorm > OptimizeConstant.RESIDUAL
       expVokBool = 0;
    end
    if (expVokBool == 0)
            LL = OptimizeConstant.LL_ERROR_VALUE;
            grad = ones(Op.n,1);
            disp('The parameters not fesible')
            return; 
    end
    % Get gradient
    gradExpV = objArray(Op.n);
    for i = 1:Op.n
        u = M .* (AttLc(i).Value); 
        v = sparse(u * Z); 
        p = Atts(i).Value(:,lastIndexNetworkState+1 : maxDest) .* Mfull(:,lastIndexNetworkState+1 : maxDest);
        p(lastIndexNetworkState+1,:) = sparse(zeros(1,maxDest - lastIndexNetworkState));        
        p = sparse(p);
        gradExpV(i).value =  sparse(A\(v + p)); 
    end
    
    % Compute the LL and gradient.
    
    for n = 1:nbobs    
%        n
        dest = Obs(n, 1);
        orig = Obs(n, 2);
        expV = Z(:,dest - lastIndexNetworkState);
        expV = full(abs(expV));         
        lnPn = - 1 * (1/mu) * log(expV(orig));
        for i = 1: Op.n
            Gradient(n,i) = - gradExpV(i).value(orig,dest - lastIndexNetworkState)/expV(orig);
        end
        
        sumInstU = 0;
        sumInstX = zeros(1,Op.n);
        
        path = Obs(n,:);
        lpath = size(find(path),2);
        % Compute regular attributes
        for i = 2:lpath - 1
            sumInstU = sumInstU + Ufull(path(i),path(i+1));
            for j = 1:Op.n
                sumInstX(j) = sumInstX(j) + Atts(j).Value(path(i),path(i+1));
            end
        end
        Gradient(n,:) = Gradient(n,:) + sumInstX;
        lnPn = lnPn + (1/mu)*sumInstU ;  
        LL =  LL + (lnPn - LL)/n;
        grad = grad + (Gradient(n,:) - grad)/n;
        Gradient(n,:) = - Gradient(n,:);
    end
    LL = -1 * LL; % IN ORDER TO HAVE A MIN PROBLEM
    grad =  - grad';
end

%%
