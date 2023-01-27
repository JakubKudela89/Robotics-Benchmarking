% main file used to run all the experiments

clear;clc;close all;
addpath("algos/LSHADE");addpath("algos/DE"); addpath("algos/PSO"); addpath("algos/AGSK");   %paths to algorithms
addpath("algos/CMA-ES"); addpath("algos/ABC"); addpath("algos/HSES"); addpath("algos/RS");

load Pos_pnts.mat;                      % collection of reachable points

angles = [0, 0.0, 0.0, 0.0, 0.0, 0.0];  % starting position of the arm

nr_runs = 20;           %number of independent runs
nr_problems = 10;       %number of instances 

for nr_points = 3:6     %3:6
for nr_changes = 1:5    %1:5

res_LSHADE = cell(nr_problems, nr_runs);
res_DE =  cell(nr_problems, nr_runs);
res_PSO =  cell(nr_problems, nr_runs);
res_ABC =  cell(nr_problems, nr_runs);
res_AGSK =  cell(nr_problems, nr_runs);
res_HSES =  cell(nr_problems, nr_runs);
res_CMAES =  cell(nr_problems, nr_runs);
res_RS =  cell(nr_problems, nr_runs);

for i=1:nr_problems
    rng(i,'Twister');                          % seed for reproducibility
    pnts = Pos_pnts(randi(1e5,nr_points,1),:); % random selection of reachable points
    [m,n] = size(pnts);

    fhd = @(x) obj_f(x, angles, pnts);          % objective function reduction

    parfor j=1:nr_runs  
        problem_size = 6*nr_changes;            % problem dimension
        max_nfes = 10000*problem_size;          % max function evaluations
        pop_size = 10*problem_size;             % common population size
        optimum = 0;                            % some algorithms need a lower bound on the optimum
        lb = -2*pi; ub = 2*pi;                  % lower and upper bounds on the variables

        fprintf('runing LSHADE, problem %u, run %u, ',i,j);
        [best_val, best_sol, res_LSHADE_str] = run_lshade(fhd,problem_size,max_nfes,pop_size,optimum,lb,ub);
        fprintf('result LSHADE %e \n',best_val);
        res_LSHADE{i,j} = res_LSHADE_str;

        fprintf('runing DE, problem %u, run %u, ',i,j);
        [res_DE_str] = Run_DE(fhd,problem_size,max_nfes,pop_size,optimum,lb,ub); 
        fprintf('result DE %e \n',res_DE_str.bestval);
        res_DE{i,j} = res_DE_str;

        fprintf('runing PSO, problem %u, run %u, ',i,j);
        [res_PSO_str]= PSO_func(problem_size,pop_size,fhd,ceil(max_nfes/pop_size),lb,ub);
        fprintf('result PSO %e \n',res_PSO_str.bestval);
        res_PSO{i,j} = res_PSO_str;

        fprintf('runing ABC, problem %u, run %u, ',i,j);
        [res_ABC_str, X, FVAL, Foods,ObjVal,BestFVALCycle] = ABC(fhd,lb*ones(1,problem_size),ub*ones(1,problem_size),max_nfes,pop_size);
        fprintf('result ABC %e \n',FVAL);
        res_ABC{i,j} = res_ABC_str;


        fprintf('runing CMAES, problem %u, run %u, ',i,j);
        res_CMAES_str = cmaes(fhd_CMAES,problem_size,max_nfes,lb,ub);
        fprintf('result CMAES %e \n',res_CMAES_str.bestval);
        res_CMAES{i,j} = res_CMAES_str;

        fprintf('runing AGSK, problem %u, run %u, ',i,j);
        res_AGSK_str = AGSK(fhd,problem_size,max_nfes,lb,ub);
        fprintf('result AGSK %e \n',res_AGSK_str.bestval);
        res_AGSK{i,j} = res_AGSK_str;

        fprintf('runing HSES, problem %u, run %u, ',i,j);
        res_HSES_str = HSES(fhd,problem_size,max_nfes,lb,ub);
        fprintf('result HSES %e \n',res_HSES_str.bestval);
        res_HSES{i,j} = res_HSES_str;

        fprintf('runing RS, problem %u, run %u, ',i,j);
        res_RS_str = rs(fhd,problem_size,max_nfes,lb,ub); 
        fprintf('result RS %e \n',res_RS_str.bestval);
        res_RS{i,j} = res_RS_str;
  end
end
str = ['res_',num2str(nr_points),'_pnts_',num2str(nr_changes),'_changes','.mat'];
save(str);      % saving results into a .mat file
end
end