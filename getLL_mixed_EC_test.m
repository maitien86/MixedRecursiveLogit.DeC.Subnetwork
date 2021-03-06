% Quasi mixed recursive logit estimation
% Compute the loglikelohood value and its gradient.
%%
function [LL, grad] = getLL_mixed_EC_test()

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
    global Omega;
    global nDraws;
    global TTSub;
    %% If IRN
    
    %% Setting parameters
    [lastIndexNetworkState, maxDest] = size(incidenceFull);
    beta = zeros(Op.n,1);
    LL = 0;
    grad = zeros(1,Op.n);    
    AttLc = objArray(Op.n);
    for i = 1:Op.n
        AttLc(i).value =  (Atts(i).value(1:lastIndexNetworkState,1:lastIndexNetworkState));
        AttLc(i).value(:,lastIndexNetworkState+1) = sparse(zeros(lastIndexNetworkState,1));
        AttLc(i).value(lastIndexNetworkState+1,:) = sparse(zeros(1, lastIndexNetworkState + 1));
    end
    % Mixed logit computation
    Zod = zeros(nDraws, nbobs);
    temp = zeros(nbobs,Op.n);
    gradExpVod = objArray(nDraws); 
    gradExpV = objArray(Op.n);
    for r = 1:nDraws
        tOmega = ones(Op.n,1);
        tOmega(Op.natt+1: Op.n) = Omega(r,:)';  
        beta = Op.x .* tOmega;
        Mfull = getM(beta, isLinkSizeInclusive);
        MregularNetwork = Mfull(1:lastIndexNetworkState,1:lastIndexNetworkState);
        Ufull = getU(beta, isLinkSizeInclusive);
        UregularNetwork = Ufull(1:lastIndexNetworkState,1:lastIndexNetworkState);
        % Initialize
        M = MregularNetwork;
        M(:,lastIndexNetworkState+1) = sparse(zeros(lastIndexNetworkState,1));
        M(lastIndexNetworkState+1,:) = sparse(zeros(1, lastIndexNetworkState + 1));
        U = UregularNetwork;
        U(:,lastIndexNetworkState+1) = sparse(zeros(lastIndexNetworkState,1));
        U(lastIndexNetworkState+1,:) = sparse(zeros(1, lastIndexNetworkState + 1));
        M = sparse(M);
        U = sparse(U);
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
        K = size(B,2);
        V = zeros(N,K*Op.n);
        for i = 1:Op.n
            v = (M .* (AttLc(i).value)) * Z * tOmega(i); 
            %v = sparse(u * Z); 
            p = tOmega(i) * Atts(i).value(:,lastIndexNetworkState+1 : maxDest) .* Mfull(:,lastIndexNetworkState+1 : maxDest);
            p(lastIndexNetworkState+1,:) = sparse(zeros(1,maxDest - lastIndexNetworkState));        
            p = sparse(p);
            V(:,(i-1)*K+1:i*K) = v+p;          
        end
        X = A\V;
        % compute the gradients of the value function
        for i = 1:Op.n
            gradExpV(i).value =  sparse(X(:,(i-1)*K+1:i*K)); 
        end
        
        % Save to store :)
        for n = 1:nbobs
            dest = Obs(n, 1);
            orig = Obs(n, 2);
            Zod(r,n) = Z(orig,dest - lastIndexNetworkState);
            for i = 1: Op.n
                temp(n,i) = gradExpV(i).value(orig,dest - lastIndexNetworkState);
            end
             
        end
        gradExpVod(r).value = temp;
    end     
    % Compute the LL and gradient for mixed logit, common random number .
    Usave = objArray(nDraws);
    for r = 1: nDraws
        tOmega = ones(Op.n,1);
        tOmega(Op.natt+1: Op.n) = Omega(r,:)';  
        beta = Op.x .* tOmega;
        Usave(r).value = sparse(getU(beta, isLinkSizeInclusive));
    end
    for n = 1:nbobs
        Pn = 0;
        gradPn = zeros(1,Op.n);
        for r = 1: nDraws      
            tOmega = ones(Op.n,1);
            tOmega(Op.natt+1: Op.n) = Omega(r,:)';            
            Ufull = Usave(r).value;
            expV0 = Zod(r,n);
            expV0 = abs(expV0);
            lnPnr = - 1 * log(expV0);
            gradLnPnr = zeros(1, Op.n);% gradient of logarite of choice probability
            for i = 1: Op.n
                gradLnPnr(i) = - gradExpVod(r).value(n,i)/expV0;
            end
            sumInstU = 0;
            sumInstX = zeros(1,Op.n);
            path = Obs(n,:);
            lpath = size(find(path),2);
            % Compute regular attributes
            for i = 2:lpath - 1
                sumInstU = sumInstU + Ufull(path(i),path(i+1));
                for j = 1:Op.n
                    sumInstX(j) = sumInstX(j) +  tOmega(j) * Atts(j).value(path(i),path(i+1));
                end
            end
            lnPnr = lnPnr + sumInstU ; 
            gradLnPnr = gradLnPnr + sumInstX;
            gradPnr = gradLnPnr * exp(lnPnr);
            Pn = Pn + (exp(lnPnr) - Pn)/r; 
            gradPn = gradPn + (gradPnr - gradPn)/r;
        end
        LL =  LL + (log(Pn) - LL)/n;
        Gradient(n,:) = gradPn / Pn;
        grad = grad + (Gradient(n,:) - grad)/n;
        Gradient(n,:) = - Gradient(n,:);
    end
    LL = -1 * LL; % IN ORDER TO HAVE A MIN PROBLEM
    grad =  - grad';
end

%%
