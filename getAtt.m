%   Get attribute 
function getAtt()

    global incidenceFull;
    global EstimatedTime;
    global TurnAngles;
    global LeftTurn;
    global Uturn;
    global LSatt;
    global isLinkSizeInclusive;
    global Atts;
    global Op;
    global TTSub;
    global corrType;
    mu = 1;
    Incidence = 1/mu * incidenceFull;
    Atts = objArray(Op.n);
    Atts(1).value = sparse(Incidence .* EstimatedTime);
    Atts(2).value = sparse(Incidence .* TurnAngles);
    Atts(3).value = sparse(Incidence .* LeftTurn);
    Atts(4).value = sparse(Incidence .* Uturn);
    if isLinkSizeInclusive == 1
        Atts(5).value = sparse(Incidence .* LSatt);
    end
    
    if  corrType == OptimizeConstant.EC
        for i = Op.natt+1 : Op.n
            Atts(i).value = TTSub(i - Op.natt).value;
        end
    end
end
