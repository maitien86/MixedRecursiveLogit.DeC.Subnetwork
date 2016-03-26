%% select hessian for retrospective model
function ind = swretro_select_hessian()
    global Op;
    ind = 0;
    mError = bitmax;
    for i =1:Op.m
        p = Op.prev_x - Op.x; 
       error = abs(Op.prev_value - get_quadratic_model(Op.value, Op.grad, p, Op.Hi(i).value))
       if (mError >= error)
           ind = i;
           mError = error;
       end
    end
    if ind ~= Op.prevHind 
        Op.nSwitch = Op.nSwitch + 1;
    end    
    Op.prevHind = ind;
    
end
%% 
%% Compute: f + g'p + 0.5p'Hp
function model = get_quadratic_model(f,g,p,H)
     model = f + g'*p + 0.5*p' * H * p;
end