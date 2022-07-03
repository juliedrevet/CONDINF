function helper_plot_corr(cfg)

if isfield(cfg,'pbar'),pbar = cfg.pbar;
else
    pbar = 1;
end
if isfield(cfg,'HGT'),HGT = cfg.HGT;
else
    HGT = 4;
end
if isfield(cfg,'mkr_siz'),mkr_siz = cfg.mkr_siz;
else
    mkr_siz = 4;
end
if isfield(cfg,'figname')
    figure('Color','white','Name',cfg.figname);
else    
    figure('Color','white');
end
if ~isfield(cfg,'rgb')
    cfg.rgb = [0.5 0.5 0.5];
end
if ~isfield(cfg,'xlabel')
    cfg.xlabel = 'data1';
end
if ~isfield(cfg,'ylabel')
    cfg.ylabel = 'data2';
end
if isfield(cfg,'verbose'),verbose = cfg.verbose;
else
    verbose = true;
end
if isfield(cfg,'identity'),identity = cfg.identity;
else
    identity = true;
end   
if isfield(cfg,'ctype'),ctype = cfg.ctype;
else
    ctype = 'Pearson';
end  

hold on

if isfield(cfg,'xlim'),xlim(cfg.xlim);end
if isfield(cfg,'ylim'),ylim(cfg.ylim);end


if isfield(cfg,'ci') && cfg.ci && isfield(cfg,'xavg') && isfield(cfg,'yavg') % add confidence interval
    rgb = cfg.rgb;
    ft = fittype('a*x + b');
    
    X = cfg.xavg;
    Y = cfg.yavg;
    [f,~] = fit(X,Y,ft,fitoptions( ...
        'Method','NonlinearLeastSquares', ...
        'StartPoint',[0,0]));
    
    xvs = linspace(min(xlim),max(xlim),100);
    yvs = feval(f,xvs)';
    yci = predint(f,xvs,0.95,'functional')';
    patch([xvs,fliplr(xvs)],[yci(1,:),fliplr(yci(2,:))],0.5*(rgb+1),'EdgeColor','none');
    plot(xvs,yvs,'-','Color',rgb,'LineWidth',2);
else
    if identity
        plot(xlim,xlim,'k')
    end
end


if isfield(cfg,'xavg') && isfield(cfg,'yavg') && isfield(cfg,'xstd') && isfield(cfg,'rgb')
    rgb = cfg.rgb;
    X = cfg.xavg;
    Y = cfg.yavg;
    ndata = numel(X);
    X_std = X + cfg.xstd*[-1,+1];
    for idata = 1:ndata
        plot(X_std(idata,:),Y(idata)*[1;1],'-','Color',(1+rgb)/2,'LineWidth',0.5);
    end
end
if isfield(cfg,'yavg') && isfield(cfg,'xavg') && isfield(cfg,'ystd') && isfield(cfg,'rgb')
    rgb = cfg.rgb;
    X = cfg.xavg; 
    Y = cfg.yavg;
    ndata = numel(Y);
    Y_std = Y + cfg.ystd*[-1,+1];
    for idata = 1:ndata
        plot(X(idata)*[1;1],Y_std(idata,:),'-','Color',(1+rgb)/2,'LineWidth',0.5);
    end
end
if isfield(cfg,'xavg')&& isfield(cfg,'yavg')
    rgb = cfg.rgb;
    X = cfg.xavg; 
    Y = cfg.yavg;
    ndata = numel(Y);
    for idata = 1:ndata
        plot(X(idata),Y(idata),'wo','MarkerFaceColor',rgb,'MarkerSize',mkr_siz,'LineWidth',0.5);
    end
end
if isfield(cfg,'xavg')&& isfield(cfg,'yavg') && verbose
    X = cfg.xavg; 
    Y = cfg.yavg;
    [rho,pval] = corr(X,Y,'type',ctype);
    fprintf('\t%s versus %s\n\t\tr = %+.3f, p = %.3f, r^2 = %+.3f\n',cfg.xlabel,cfg.ylabel,rho,pval,rho^2);
    if pval <= 0.05
    text(.05*(max(xlim)-min(xlim))+min(xlim),.9*(max(ylim)-min(ylim))+min(ylim),sprintf('r^2 = %4.3f\np < %4.2g\n',rho^2,ceil(1000*pval)/1000))
    else
        text(.05*(max(xlim)-min(xlim))+min(xlim),.9*(max(ylim)-min(ylim))+min(ylim),sprintf('r^2 = %4.3f\np > %4.2g\n',rho^2,floor(1000*pval)/1000))
    end
end

if isfield(cfg,'xtick')
    set(gca,'XTick',cfg.xtick);
end
if isfield(cfg,'xticklabel')
    set(gca,'XTickLabel',cfg.xticklabel);
end
 if isfield(cfg,'ytick')
    set(gca,'YTick',cfg.ytick);
end
if isfield(cfg,'yticklabel')
    set(gca,'YTickLabel',cfg.yticklabel);
end   

set(gca,'Layer','top','Box','off','TickDir','out','PlotBoxAspectRatio',[pbar,1,1]);
xlabel(cfg.xlabel);
ylabel(cfg.ylabel);

fig = helper_ILL(gcf,pbar,HGT,mkr_siz);
set(fig,'PaperPositionMode','manual', ...
    'PaperPosition',[2.5,13,16,HGT],'PaperUnits','centimeters', ...
    'PaperType','A4','PaperOrientation','portrait');
if isfield(cfg,'saveornot') && isfield(cfg,'figname') && isfield(cfg,'figpath')
    if cfg.saveornot
        print(fullfile(cfg.figpath,cfg.figname),'-painters','-dpdf');
        fprintf('\nfigure %s saved in %s !\n\n',cfg.figname,cfg.figpath)
    end
end
        
end