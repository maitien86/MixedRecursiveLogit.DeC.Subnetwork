%   Multi Hessians approximation
%   BHHH - BFGS - SR1 - BIGGS - Compbined_BFGS - Conbined_SR1
%%
function ok = update_multi_hessians_approx()
    global Op;
    for i = 1:Op.m
        switch Op.set_hes(i)
            case OptimizeConstant.BHHH
                Op.Hi(i).value = BHHH();
                ok = true;
            case OptimizeConstant.SR1
                [Op.Hi(i).value, ok] = SR1(Op.step, Op.deltaGrad, Op.Hi(i).value);
            case OptimizeConstant.BFGS
                [Op.Hi(i).value, ok] = BFGS(Op.step, Op.deltaGrad, Op.Hi(i).value);            
            case OptimizeConstant.CB_BFGS
                [Op.Hi(i).value, ok] = Combined_BFGS(Op.step, Op.deltaGrad, Op.Ai(i).value);
            case OptimizeConstant.CB_SR1
                [Op.Hi(i).value, ok] = Combined_SR1(Op.step, Op.deltaGrad, Op.Ai(i).value);
            case OptimizeConstant.BIGGS
                [Op.Hi(i).value, ok] = BIGGS(Op.deltaValue, Op.step, Op.deltaGrad, Op.grad, Op.Hi(i).value);
            case OptimizeConstant.SSA_BFGS
                [Op.Hi(i).value, ok] = SecantStatistical_approx(Op.step, Op.deltaGrad, OptimizeConstant.SSA_BFGS);
            case OptimizeConstant.SSA_SR1
                [Op.Hi(i).value, ok] = SecantStatistical_approx(Op.step, Op.deltaGrad, OptimizeConstant.SSA_SR1);
        end
    end
end