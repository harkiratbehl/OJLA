function [x,v] = gurobi_ilp(f, A, b, Aeq, beq)
try
    clear model;
    model.obj = double(f);%costofalledges
    model.A = sparse([A;Aeq]);%(bb+trgts)*edges
    model.rhs = [b;beq];%(bb+trgts)*1
    [Ah,~] = size(A);%noofbb
    [Aeqh,~] = size(Aeq);%nooftrgts
    model.sense = char(['<' * ones(1, Ah), '=' * ones(1,Aeqh)]);
    model.vtype = 'B';
    model.modelsense = 'min';
    clear params;
    params.outputflag = 0;
    result = gurobi(model, params);
    x = result.x;
    v = result.objval;
catch gurobiError
    x = [];
    v = [];
end
end

