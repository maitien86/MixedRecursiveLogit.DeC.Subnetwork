% Quasi mixed recursive logit estimation
% Compute the loglikelohood value and its gradient.
%%
function [LL, grad] = getLL_mixed_IRN()

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
    % ----------------------------------------------------
    % Setting parameters
    [lastIndexNetworkState, maxDest] = size(incidenceFull);
    beta = Op.x(2:end,1);
    LL = 0;
    grad = zeros(1,Op.n);    
    AttLc = objArray(Op.n);
    for i = 2:Op.n
        AttLc(i).value =  (Atts(i-1).Value(1:lastIndexNetworkState,1:lastIndexNetworkState));
        AttLc(i).value(:,lastIndexNetworkState+1) = sparse(zeros(lastIndexNetworkState,1));
        AttLc(i).value(lastIndexNetworkState+1,:) = sparse(zeros(1, lastIndexNetworkState + 1));
    end
    AttLc(1).value = AttLc(2).value;
    % Mixed logit computation
    gradExpV = objArray(Op.n);
    for n = 1: nbobs
        dest = Obs(n, 1);
        orig = Obs(n, 2);
        Pn = 0;
        gradPn = zeros(1,Op.n);
        for r = 1:nDraws
            beta(1) = Op.x(1) + Op.x(2) * Omega(n, r);
            AttLc(2).value = AttLc(1).value * Omega(n, r);    
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
            % Compute matrix Z
            N = size(M,1);
            B = sparse(zeros(N,1));
            B(N) = 1;
            B(1:lastIndexNetworkState,1) = Mfull(:,dest);
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
            for i = 1:Op.n
                u = M .* (AttLc(i).value); 
                v = sparse(u * Z); 
                if i > 2
                    p = Atts(i-1).Value(:,dest) .* Mfull(:,dest);
                elseif i == 2
                    p = Omega(n, r) * Atts(i-1).Value(:,dest) .* Mfull(:,dest);
                else
                    p = Atts(1).Value(:,dest) .* Mfull(:,dest);  
                end
                p(lastIndexNetworkState+1) = 0;        
                p = sparse(p);
                gradExpV(i).value =  sparse(A\(v + p)); 
            end

            
            %Ufull = Usave(r).value;
            expV0 = Z(orig);
            expV0 = abs(expV0);
            lnPnr = - 1 * log(expV0);
            gradLnPnr = zeros(1, Op.n);% gradient of logarite of choice probability
            for i = 1: Op.n
                gradLnPnr(i) = - gradExpV(i).value(orig)/expV0;
            end
            sumInstU = 0;
            sumInstX = zeros(1,Op.n);
            path = Obs(n,:);
            lpath = size(find(path),2);
            % Compute regular attributes
            for i = 2:lpath - 1
                sumInstU = sumInstU + Ufull(path(i),path(i+1));
                for j = 3:Op.n
                    sumInstX(j) = sumInstX(j) + Atts(j-1).Value(path(i),path(i+1));
                end
                sumInstX(2) = sumInstX(2) + Atts(1).Value(path(i),path(i+1)) * Omega(n, r);
                sumInstX(1) = sumInstX(1) + Atts(1).Value(path(i),path(i+1));
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
