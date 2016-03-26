% Initialie optimization structure for predictive model
%%
%
function [] = initialize_switching_structure()
    global Op;
    Op.Switching = OptimizeConstant.SWITCHING;
    if strcmp(Op.Optim_Method,OptimizeConstant.TRUST_REGION_METHOD)
        Op.m = 3;
    else
        Op.m = 2;
    end
    Op.set_hes = [OptimizeConstant.BHHH;OptimizeConstant.CB_BFGS; OptimizeConstant.CB_SR1;  OptimizeConstant.SSA_BFGS; OptimizeConstant.SSA_SR1; ];
    %Op.set_hes = [OptimizeConstant.BHHH;OptimizeConstant.CB_BFGS];
    Op.Hi = objArray(Op.m);
    Op.Ai = objArray(Op.m);
    for i = 1:Op.m
        Op.Hi(i).value = eye(Op.n);
        if Op.set_hes(i) == OptimizeConstant.CB_SR1 || Op.set_hes(i) == OptimizeConstant.CB_BFGS
            Op.Ai(i).value = zeros(Op.n);
        end
    end
end