function [mean_act] = receptor_output(env,rvec,param)
%receptor activity given environment, circuit parameters, and receptor
%arrangement

arguments
    env
    rvec
    param
end

if param.noisy
    c_samp = poissrnd(repmat(env,param.nsamp,1));
    a_samp = poissrnd(rvec.*hillfun(c_samp,param)/param.rtot);
    mean_act = mean(a_samp,1);
else
    mean_act = rvec.*hillfun(env,param)/param.rtot;
end

end

