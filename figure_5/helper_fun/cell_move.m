function [results] = cell_move(param,varargin)

p = inputParser;
addRequired(p,'param',@isstruct);
addParameter(p,'verbose',false,@islogical);
addParameter(p,'plotdynamics',false,@islogical);
addParameter(p,'mode',"static",@isstring);
% options: static, dynamic, localization (if static, decoder_method not
% needed)
addParameter(p,'decoder_method',"optimal_noise",@isstring);
% options: optimal_noise, perfect, randomwalk, yaxis1d
addParameter(p,'hasMemory',false,@islogical);
% whether to allow cell to perform temporal averaging
addParameter(p,'receptor',"uniform",@isstring);
% options: feedback, uniform, w1dist

parse(p,param,varargin{:});
param = p.Results.param;
isVerbose = p.Results.verbose;
plotdynamics = p.Results.plotdynamics;
mode = p.Results.mode;
decoder_method = p.Results.decoder_method;
hasMemory = p.Results.hasMemory;
receptor = p.Results.receptor;

it = 1;
m = param.N;

cellp = param.cellp;
stepsz = param.stepsz;
fcount = param.fcount;
time_size = floor(param.T/param.dt);
move_rate = floor(30/param.dt); %move cell every 30 seconds
s = param.mean_cell_radius;

% Preallocate space
results.f = zeros(time_size,m);
results.cellp = zeros(time_size,2);
results.env = zeros(time_size,m);
if hasMemory
    memlen = move_rate*10; % [seconds] memory length
    results.an = zeros(memlen,m);
end

% Initial environment
results.cellp(it,:) = cellp;
cellboundary = param.cellboundary;
coord = cellp + cellboundary;
env = arrayfun(fcount,coord(:,1),coord(:,2))';
results.env(it,:) = env;
format long

if isequal(receptor, 'feedback')
    [param.L, param.R] = CNMatrix(param);
    results.FC = zeros(time_size,1);
    results.stat = zeros(1,4); 
    results.kfb = zeros(time_size,1);
    results.f(it,:)  = 0.9*param.rtot/m*ones(1,m);
    results.FC(it,:) = param.rtot - sum(results.f(it,:));
    %initialize param.hprop
    rvec = results.f(it,:);
    param.hprop = param.h/(mean(receptor_output(env,rvec,param)));
    results.kfb(it) = param.h;
    
elseif isequal(receptor, 'uniform')
    results.f(it,:) = param.rtot/m * ones(1,m);
    
elseif isequal(receptor,'w1dist')
    epsilon = 0.1/m;
    results.f(it,:) = (-epsilon + (epsilon-(-epsilon)).*rand(1,m) + 1/m)*param.rtot;
    gamma = param.gamma;
    % Optimization parameter
    problem.options = optimoptions('fmincon','Algorithm','sqp',...
        'SpecifyObjectiveGradient',true,...
        'Display','off','CheckGradient',false);
    problem.solver = 'fmincon';
    problem.Aeq = ones(1,m);
    problem.beq = 1;
    problem.lb = zeros(1,m);
    problem.ub = ones(1,m);
    problem.x0 = results.f(it,:)./sum(results.f(it,:));
end

for it = 2:time_size  %it == "time iterators"
    % update receptors
    if isequal(receptor,'feedback')
        results = update_receptor(param,env,results,it);
    elseif isequal(receptor,'uniform')
        results.f(it,:) = results.f(1,:);
    elseif isequal(receptor,'w1dist')
        if mod(it,move_rate) == 1
            exlogx = compute_exlogx(env,param);
            problem.objective = @(x) totalMI(x,env,exlogx,param,...
                'gamma',gamma,'lastrvec',results.f(it-move_rate,:));
            [xopt,~] = fmincon(problem);
            results.f(it,:) = xopt*param.rtot;
        else
            results.f(it,:) = results.f(it-1,:);
        end
    end
    if hasMemory
        an = receptor_output(env,results.f(it-1,:),param);
        results.an = [an; results.an(1:min(it,memlen)-1,:)];
    end
    
    % move cell and update environment
    if ~isequal(mode,'static')
        if mod(it,move_rate) == 1
            rvec = results.f(it,:);
            ansum = receptor_output(env,rvec,param);
            if hasMemory
                ansum = mean([ansum; results.an(1:min(it,memlen)-1,:)]);
            end
            movdir = grad_decode(ansum, decoder_method);
            
            %check if cell still inbound
            cellp2 = cellp + stepsz*[cos(movdir), sin(movdir)];
            if inbound(cellp2 + cellboundary) || isequal(decoder_method,"yaxis1d")
                cellp = cellp2; %change cell position
                coord = cellp + cellboundary;
                env = arrayfun(fcount,coord(:,1),coord(:,2))';
            end
            if isequal(mode,'localization') 
                if cellp(1) < 5+s % stop when cell is within 5um of source
                    results.cellp(it,:) = cellp;
                    results.env(it,:) = env;
                    if isequal(receptor,'feedback') %record summary statistics
                        results.stat = record_stat(param,results,it);
                    end
                    if isVerbose
                        disp(strcat('(',receptor,') Time taken: ', num2str(floor(it*param.dt/60)),' mins'))
                    end
                    break
                end
            end
            if plotdynamics
                plotting_dynamics(rvec,env,it,param,cellp);
            end
        end
        
    end
    
    
    results.cellp(it,:) = cellp;
    results.env(it,:) = env;
    
    if it == time_size %record summary statistics
        if isequal(receptor,'feedback')
            results.stat = record_stat(param,results,it);
        end
        if isequal(mode,'localization') && isVerbose
            disp(strcat("(",receptor,") ", "unfinished"))
        end
    end
end

end

function plotting_dynamics(rvec,env,it,param,cellp)
figure(3)
yyaxis left
plot(rvec,'linewidth',2)
ylabel('Receptor distribution','fontsize',16)
yyaxis right
plot(env,'linewidth',2)
ylabel('Ligand concentration','fontsize',16)
title(strcat("Time elapsed = ", num2str(floor(it*param.dt/60))," mins",...
       " cell pos = ", num2str(cellp(1)),', ',num2str(cellp(2))),'fontsize',14);
end
