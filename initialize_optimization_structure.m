%   Initialize optimization structure
%%
function [] = initialize_optimization_structure()
    global Op;
    global isLinkSizeInclusive;
    global isFixedUturn;
    global nbobs;
    global corrType;
    global nECs;
    Op.Optim_Method = OptimizeConstant.TRUST_REGION_METHOD;
    Op.ETA1 = 0.05;
    Op.ETA2 = 0.75;
    Op.maxIter = 300;
    Op.k = 0;
    if corrType == OptimizeConstant.RP
        Op.n = 5;
    else
        Op.n = 4 + nECs;
    end
    if isLinkSizeInclusive == true
        Op.n = Op.n + 1;
        Op.natt = 5;
    else 
        Op.natt = 4;
    end
    
    if isFixedUturn == true
        Op.n = Op.n - 1;
    end
    Op.x = -ones(Op.n,1) * 1;
    Op.prev_x = 0 * Op.x;
    Op.prev_value = 0;
    Op.tol = 1e-5;
    Op.radius = 1.0;
    Op.Ak = zeros(Op.n);
    Op.H = eye(Op.n);
    %Op.Hi = objArray(nbobs);
    %for i = 1:nbobs
    %    Op.Hi(i).value = eye(Op.n);
    %end
end