function parameter_scan(task,savefname)
%% parameter phase diagram (SI)
p = inputParser;
addRequired(p,'task',@isstring);
addRequired(p,'savefname',@isstring);

parse(p,task,savefname);
task = p.Results.task; % retention or localization;
savefname = p.Results.savefname;

% setting scheme parameter
nval = 15;
kofflist = logspace(-2,-0.5,nval); 
hlist = logspace(-3.5,-2,nval); 
paramlist = combvec(kofflist,hlist)';
nparam = size(paramlist,1);

%% setting paramaters
load('tissue_300by900_szopt','receptor_params')
load('scheme_parameter')

%% localization
% setting feedback scheme parameter
if isequal(task, "localization")
    thour = 1;
    conversion_factor = conc2count(scheme_param.mean_cell_radius,scheme_param.N);
    scheme_param.T = 60*60*thour; % total time in seconds
    scheme_param.kd = receptor_params.kd*conversion_factor;
    scheme_param.rtot = receptor_params.rtot;
    scheme_param.receptornoise = receptor_params.receptornoise;
    
    scheme_succ = zeros(nparam,1);
    scheme_stat = zeros(nparam,4);
    fname = "tissue_300by900";
    for ii = 1:nparam
        param = scheme_param;
        param.koff = paramlist(ii,1);
        param.h = paramlist(ii,2)/2;
        [res, ~,stat] = racing_cells(fname, param,"envmodel","tissue",...
        "task","localization",...
        "rununif",false,...
        "save_data",0,"index",strcat("_",string(ii)));
        scheme_succ(ii) = res;
        scheme_stat(ii,:) = stat;
        if mod(ii,5) == 0
            disp(ii)
        end
    end
end

%% retention
if isequal(task, 'retention')
    thour = 1;
    scheme_param.mean_cell_radius = ellipsePerimeter(2,5)/(2*pi);
    conversion_factor = conc2count(scheme_param.mean_cell_radius,scheme_param.N);
    scheme_param.T = 60*60*thour; % seconds
    scheme_param.kd = receptor_params.kd*conversion_factor;
    scheme_param.rtot = receptor_params.rtot;
    scheme_param.receptornoise = receptor_params.receptornoise;
    
    scheme_err = zeros(nparam,1);
    scheme_stat = zeros(nparam,4);
    fname = "tissue_retention_1";
    for ii = 1:nparam
        param = scheme_param;
        param.koff = paramlist(ii,1);
        param.h = paramlist(ii,2);
        [res, ~, stat] = racing_cells(fname, param,"envmodel","tissue",...
            "save_data",0,...
            "rununif",false,...
            "task","retention","index",strcat("_",string(ii)));
        scheme_err(ii) = res;
        scheme_stat(ii,:) = stat;
        if mod(ii,5) == 0
            disp(ii)
        end
    end
end

save(savefname);

end
