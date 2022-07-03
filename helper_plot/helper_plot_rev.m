function helper_plot_rev(cfg)
% this funciton plots reversal curves
% possible fields:
%       * yavg: plain lines, dim1 = ntrials around REV,dim2 = nmodel to plot
%       * yerr: shaded areas sem, dim1 = ntrials around REV,dim2 = nmodel to plot
%       * yavg_dashed: dashed gray lines, dim1 = ntrials around REV,dim2 = nmodel to plot
%       * yerr_dashed: shaded gray areas sem, dim1 = ntrials around REV,dim2 = nmodel to plot
%       * Yavg : data (circles)

if isfield(cfg,'pbar'),pbar = cfg.pbar;
else
    pbar = 1.1;
end
if isfield(cfg,'HGT'),HGT = cfg.HGT;
else
    HGT = 4;
end
if isfield(cfg,'mkr_siz'),mkr_siz = cfg.mkr_siz;
else
    mkr_siz = 4;
end

% models to plot in color, plain lines
if isfield(cfg,'yavg')
    yavg = cfg.yavg;
    nmodel = size(yavg,2);
end
if isfield(cfg,'rgb'),rgb = cfg.rgb;
else, rgb = lines; rgb = rgb(1:nmodel,:);
end
if isfield(cfg,'yavg')&&(nmodel>size(rgb,1))
    rgb = cat(1,rgb,lines);
end


% models to plot in gray, dashed lines
if isfield(cfg,'yavg_dashed')
    yavg_dashed = cfg.yavg_dashed;
    if any(size(yavg_dashed)==1) % only one model, transform to column
        yavg_dashed = yavg_dashed(:);
    end
    nmodel_dashed = size(yavg_dashed,2);
    if isfield(cfg,'rgb_dashed'),rgb_dashed = cfg.rgb_dashed;
    else
        rgb_dashed = (0.8:(0.1/nmodel_dashed):0.9)'.*[1 1 1];
    end
end
if isfield(cfg,'yavg_dashed')&&(nmodel_dashed>size(rgb_dashed,1))
    rgb_dashed = repmat(rgb_dashed,1,ceil(nmodel_dashed/size(rgb_dashed,1)));
end


% participants to plot in color, plain dots
if isfield(cfg,'Yavg')
    Yavg = cfg.Yavg;
    nmodel_dat = size(Yavg,2);
end
if isfield(cfg,'rgb_dat'),rgb_dat = cfg.rgb_dat;
elseif isfield(cfg,'Yavg'), rgb_dat = lines; rgb_dat = rgb_dat(1:nmodel_dat,:);
end
if isfield(cfg,'Yavg')&&(nmodel_dat>size(rgb_dat,1))
    rgb_dat = cat(1,rgb_dat,lines);
end

if isfield(cfg,'figname')
    figname = cfg.figname;
    figure('Color','white','Name',figname);
else    
    figure('Color','white');
end

hold on

if isfield(cfg,'xvec'),xvec = cfg.xvec;
else
    xvec = -3.5:+3.5;
end

if isfield(cfg,'xlim'),xlim(cfg.xlim);
else
    xlim([-4,+4]);
end

if isfield(cfg,'ylim'),ylim(cfg.ylim);
else
    ylim([0,1]);
end

%% SHADED AREAS
% dashed lines
if isfield(cfg,'yavg_dashed') && isfield(cfg,'yerr_dashed')
    yerr_dashed = cfg.yerr_dashed;
    nmod_err_dashed = size(yerr_dashed,2);
    nmod = max(nmodel_dashed,nmod_err_dashed);
    ntrl = size(yavg_dashed,1);
    yerr_vec = nan(2,ntrl,nmod);
    for imod = 1:nmod
        yerr_vec(1,:,imod) = reshape(yavg_dashed(:,imod),1,ntrl) + reshape(yerr_dashed(:,imod),1,ntrl);
        yerr_vec(2,:,imod) = reshape(yavg_dashed(:,imod),1,ntrl) - reshape(yerr_dashed(:,imod),1,ntrl);
        patch([xvec,fliplr(xvec)],[yerr_vec(1,:,imod),fliplr(yerr_vec(2,:,imod))], ...
            0.5*(rgb_dashed(imod,:)+1),'EdgeColor','none');
    end
end
% model
if isfield(cfg,'yavg') && isfield(cfg,'yerr')
    yerr = cfg.yerr;
    nmod_err = size(yerr,2);
    nmod = max(nmodel,nmod_err);
    ntrl = size(yavg,1);
    yerr_vec = nan(2,ntrl,nmod);
    for imod = 1:nmod
        yerr_vec(1,:,imod) = reshape(yavg(:,imod),1,ntrl) + reshape(yerr(:,imod),1,ntrl);
        yerr_vec(2,:,imod) = reshape(yavg(:,imod),1,ntrl) - reshape(yerr(:,imod),1,ntrl);
        patch([xvec,fliplr(xvec)],[yerr_vec(1,:,imod),fliplr(yerr_vec(2,:,imod))], ...
            0.5*(rgb(imod,:)+1),'EdgeColor','none');
    end
end

%% black lines
plot([0,0],ylim,'k');
plot(xlim,0.5*[1,1],'k');
    
%% model simulations
% dashed model simulation
if isfield(cfg,'yavg_dashed')
    for imod = 1:nmodel_dashed
        plot(xvec,yavg_dashed(:,imod),'Color',rgb_dashed(imod,:),'LineWidth',1.5,'LineStyle',':');
    end
end
% stabilized model simulation
if isfield(cfg,'yavg')
    for imod = 1:nmodel
        plot(xvec,yavg(:,imod),'Color',rgb(imod,:),'LineWidth',1.5);
    end
end

%% participants 
% error bars
if isfield(cfg,'Yavg') && isfield(cfg,'Yerr')
    if all(size(cfg.Yavg)==size(cfg.Yerr))
        nmodel_dat =  size(cfg.Yavg,2);
        for imod = 1:nmodel_dat
            for i = 1:8
                plot(xvec(i)*[1,1],cfg.Yavg(i,imod)+cfg.Yerr(i,imod)*[+1,-1],'k-');
            end
        end
    end
end

% participants dots
if isfield(cfg,'Yavg')
    nmodel_dat =  size(cfg.Yavg,2);
    for imod = 1:nmodel_dat
        plot(xvec,cfg.Yavg(:,imod),'ko','MarkerSize',mkr_siz,'MarkerFaceColor',rgb_dat(imod,:));
    end
end


%%
set(gca,'Layer','top','Box','off','TickDir','out','PlotBoxAspectRatio',[pbar,1,1]);
set(gca,'XTick',xvec,'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'},'YTick',0:0.2:1);
xlabel('trial position from reversal');
ylabel('fraction reversed');
helper_ILL(gcf,pbar,HGT,4.5);

if isfield(cfg,'saveornot') && isfield(cfg,'figname') && isfield(cfg,'figpath')
    if cfg.saveornot
        print(fullfile(cfg.figpath,cfg.figname),'-painters','-dpdf');
        fprintf('\nfigure %s saved in %s !\n',cfg.figname,cfg.figpath)
    end
end
        