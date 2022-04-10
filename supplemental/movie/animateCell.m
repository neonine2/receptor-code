clear all;
close all;

%animate cell
load(which("tissue_300by900_localization_feedback.mat"),...
                'param', 'recC', 'recF', 'posScheme', 'posUnif','m',...
                'move_rate')

%% panel F snapshots
ylimit = -70:70;
xlimit = -49:19;
fconc = param.fcount;
s = param.mean_cell_radius;
nstep = 240;
nrun = size(recC,3);
position = cat(4,posScheme,posUnif);

fname = "haptotaxis";
if isfile(strcat(fname,".mp4"))
    fname = strcat(fname,"-1",".mp4");
end
vidfile = VideoWriter(fname,'MPEG-4');
vidfile.FrameRate = 20;
vidfile.Quality = 100;
open(vidfile);
load('TissueColormap')

disp(sum(sum(position(121,:,1:5:100,1))==0))
disp(sum(sum(position(121,:,1:5:100,2))==0))
for runind = 1:5:10 %[126,136,137,141,148,156,159]
    flag = 0;
    initp = posUnif(1,:,runind);
    coord = combvec(xlimit,ylimit)' + initp;
    envtissue = reshape(arrayfun(fconc,coord(:,1),coord(:,2)),...
        length(xlimit),length(ylimit))';
    tiledlayout(1,2,'TileSpacing','compact','padding','compact');
    a=nexttile;
    b=nexttile;
    tiles = [a,b];
    set(gcf, 'Position',  [100, 100, 14*(max(xlimit)-min(xlimit)),...
                                     6*(max(ylimit)-min(ylimit))])
    if sum(position(121,:,runind,1))==0
        schemeStatus = 'Completed < 1hr';
        stopScheme = find(sum(position(:,:,runind,1),2)==0,1);
    else
        schemeStatus = 'Failed';
        stopScheme = nstep;
    end
    if sum(position(121,:,runind,2))==0
        unifStatus = 'Completed < 1hr';
        stopUnif = find(sum(position(:,:,runind,2),2)==0,1);
    else
        unifStatus = 'Failed';
        stopUnif = nstep;
    end
    for ii = 1:nstep
        for jj = 1:2
            axes(tiles(jj))
            colormap(TissueColormap);
            imagesc(envtissue)
            colorbar()
            cellcenter = position(ii,:,runind,jj);
            if isequal(cellcenter,[0,0]) && isequal(unifStatus,'Fail')
                flag=1;
                break
            end
            if isequal(position(ii,:,runind,2),[0,0])
                flag=1;
                break
            end
            if(jj == 1)
                status = schemeStatus;
                receptor = recF(ii,:,runind)';
            else
                status = unifStatus;
                receptor = mean(recF(10,:,runind));
            end
            angle = linspace(0,2*pi*(1-1/m),m);
            cellcenter = cellcenter - initp - [min(xlimit),min(ylimit)];
            cellsurf = cellcenter+(s+receptor/10).*[cos(angle)',sin(angle)'];  
            circle = cellcenter + s.*[cos(angle)',sin(angle)']; 

            %plotting
            hold on
            fill(cellsurf(:,1),cellsurf(:,2),[17 17 17]/18,'Linewidth',1)
            fill(circle(:,1),circle(:,2),[17 17 17]/18,'Linewidth',1,...
                'Facecolor','#0072BD')
            hold off
            set(gca,'fontsize',14)
            if jj==1
                title("Feedback scheme",'fontsize',20)
            else
                title("Uniform receptors",'fontsize',20)
            end
            subtitle(status,'fontsize',18)
            if jj==1 && ii>stopScheme
                timestr = append(num2str(move_rate*stopScheme/60,'%3.2f')," min");
            else
                timestr = append(num2str(move_rate*ii/60,'%3.2f')," min");
            end
            text(45,8,timestr,'FontSize',23,'Color','white');
            pbaspect([max(xlimit)-min(xlimit),max(ylimit)-min(ylimit),1])
%             set(gcf,'position',[10,10,492*2,276*2])
        end
        frame = getframe(gcf);
        writeVideo(vidfile, frame);
        if(flag==1)
            break
        end
    end
end

movie2gif(vidfile, 'sin.gif', 'LoopCount', 0, 'DelayTime', 0)
close(vidfile)