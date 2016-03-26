%   Compute Link Size attribute from data
%   
%%
function ExpV_is_ok = getLinkSizeAtt()
    global incidenceFull; 
    global Obs;     % Observation
    global LSatt;
    global EstimatedTime;
    global TurnAngles;
    global LeftTurn;
    global Uturn;    
    beta = [-2.5,-1,-0.4,-4]';
    % ----------------------------------------------------
    mu = 1; % MU IS NORMALIZED TO ONE
    Obs_temp = Obs;
    Obs_temp(find(Obs)) =  Obs(find(Obs)) + 1;
    % Save data
    IF = (incidenceFull);    
    IFext = transferMarix(incidenceFull);
    EText = transferMarix(EstimatedTime);
    TAext = transferMarix(TurnAngles);
    LText = transferMarix(LeftTurn);
    UText = transferMarix(Uturn);    
    lastIndexNetworkState = size(IFext,1);
    dest = unique(Obs_temp(:,1));
    orig = unique(Obs_temp(:,2));
    IFext(1,orig') = 1;
    EText(1,orig') = 0.05;
    IFext(dest,lastIndexNetworkState) = 1;
    EText(dest,lastIndexNetworkState) = 0.05; 
    
    u =  beta(1) * EText + beta(2) * TAext +  beta(3) * LText + beta(4) * UText;
    
    expM = sparse(u);
    expM(find(IFext)) = exp(u(find(IFext)));
    M = expM .* IFext;    
    
    [expV, expVokBool] = getExpV(M);
    if (expVokBool == 0)
       ExpV_is_ok = false;
       disp('ExpV is not fesible')
       return; 
    end  
    P = getP(expV, M);
    G = sparse(zeros(size(expV)));
    G(1) = 1;
    I = speye(size(P));
    F = (I-P')\G;                        
    if (min(F) < 0)                
        ToZero = find(F <= 0);
        for i=1:size(ToZero,1)
            F(ToZero(i)) = 1e-9;
        end
    end
    ips = F;
    ips = ips(1:lastIndexNetworkState); 
    ips(size(IFext,2),1) = 0;   
    I = find(IFext);
    [nbnonzero, c] = size(I);
    ind1 = zeros(nbnonzero,1);
    ind2 = zeros(nbnonzero,1);
    s = zeros(nbnonzero,1);
    for i = 1:nbnonzero
        [k a] = ind2sub(size(IFext), I(i));
        ind1(i) = k;
        ind2(i) = a;
        s(i) =  ips(a);
    end    
    LSext = sparse(ind1, ind2, s, size(IFext,1), size(IFext,2));
    LSatt = LSext(2:size(IF,1)+1,2:size(IF,2)+1);
    ExpV_is_ok = true;
end