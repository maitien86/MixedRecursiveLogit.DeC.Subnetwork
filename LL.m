function [f g] = LL(x)
    global Op;
    global corrType;
    Op.x = x;
    if corrType == OptimizeConstant.RP
        [f,g] = getLL_mixed();
    else
        [f,g] = getLL_mixed_EC();
    end
    Op.nFev  = Op.nFev + 1;
end