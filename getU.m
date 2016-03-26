%   Get Utility
%%
function Ufull = getU(x, isLS)

    global incidenceFull;
    global EstimatedTime;
    global TurnAngles;
    global LeftTurn;
    global Uturn;  
    global LSatt;
    global Op;
    global TTSub;
    global corrType;
    
    u1 = x(1) * EstimatedTime;
    u2 = x(2) * TurnAngles;
    u3 = x(3) * LeftTurn;
    u4 = x(4) * Uturn;
    if isLS == true
        u5 = x(5) * LSatt;
    else
        u5 = 0 * LeftTurn;
    end
    u = sparse(u1 + u2 + u3 + u4 + u5);    
    if corrType == OptimizeConstant.EC
        for i = Op.natt + 1: Op.n;
           u = u + TTSub(i-Op.natt).value * x(i);
        end
    end
    
    u = sparse(u);    
    Ufull = incidenceFull .* u;
end


