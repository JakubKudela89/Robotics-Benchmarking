function res =rs(fhd,problem_size,max_nfes,lb,ub)
bestval = inf;
for i=1:max_nfes
    sol = lb + (ub-lb)*rand(1,problem_size);
    sol_val = fhd(sol);
    if sol_val < bestval
        bestval = sol_val;
        opt_sol = sol;
    end
end

res.bestval = bestval;
res.opt_sol = opt_sol;

end

