function [results] = update_receptor(param, env, results, it)

if it == floor(10/param.dt)
    rvec = results.f(it-1,:);
    avgA = mean(receptor_output(env,rvec,param));
    param.hprop = param.h/avgA;
end

%1) update Crank-Nicholson matricies
[L, R, meankfb] = UpdateCNMatrix(param,env,results,it);

%2) solve for next time step
v = [results.f(it-1,:)'; results.FC(it-1,:)];
C = L\R;
w = C*v;
results.f(it,:)  = w(1:end-1);
results.FC(it,:) = w(end);
results.kfb(it) = meankfb;

end

% -------------------------------------------------------------------------
function [L, R, meankfb] = UpdateCNMatrix(param,env,results,it)
% old L and R
L = param.L;
R = param.R;

% receptor activity
rvec = results.f(it-1,:);
an = receptor_output(env,rvec,param);

% setting parameters
kfb = param.hprop*an;
meankfb = mean(kfb);

% computing L and R
transportL   = - (kfb)/2;
transportR   = + (kfb)/2;
endocytosisL = (param.koff/2).*ones(1,param.N);
endocytosisR = -(param.koff./2).*ones(1,param.N);

L = L + diag([endocytosisL 0]);
L(1:end-1,end) = transportL;

R = R + diag([endocytosisR 0]);
R(1:end-1,end) = transportR;
end
