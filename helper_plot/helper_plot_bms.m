function helper_plot_bms(cfg)
% needed on the path:
%       * spm_BMS.m in SPM12 (Wellcome Trust Center for Human Neuroimaging; http://www.fil.ion.ucl.ac.uk/spm)

if isfield(cfg,'HGT'),HGT = cfg.HGT;
else, HGT = 4 ;end
if isfield(cfg,'rgb'),rgb = cfg.rgb;
else, rgb = lines; end
if isfield(cfg,'width'),wdth = cfg.width;
else, wdth = 0.7; end

if ~isfield(cfg,'elbo')
    error('no BMS possible, provide ELBO!')
end

elbo = cfg.elbo;

if size(elbo,3)>1
    error('ELBO should be a 2D-matrix: ndata x nmodel')
end

nmodel = size(elbo,2);
[pavg,pstd,pexc] = run_bms(elbo);

if isfield(cfg,'pbar'),pbar = cfg.pbar;
else, pbar = (nmodel+1.2)/HGT ;end

if isfield(cfg,'models'),models = cfg.models;
else, models = cellfun(@(s)sprintf('M%s',s),cellstr(num2str((1:nmodel)')),'UniformOutput',false);end

if isfield(cfg,'figname')
    figure('Color','white','Name',cfg.figname);
else    
    figure('Color','white');
end

hold on

if isfield(cfg,'xlim'), xlim(cfg.xlim);
else, xlim([0.4,nmodel+.6]);end
if isfield(cfg,'ylim'), ylim(cfg.ylim);
else, ylim([0 1]);end

    [~,idxmax] = max(pavg);
    plot(xlim,1/nmodel*[1,1],'k:');
    for i = 1:nmodel
        bar(i,pavg(i),wdth,'EdgeColor',rgb(i,:),'FaceColor',0.5*(rgb(i,:)+1),'LineWidth',1);
        plot(i*[1,1],pavg(i)+pstd(i)*[-1,+1],'k-');
    end
    text(idxmax,.92,sprintf('p_{exc} > %3.3f',floor(1000*pexc(idxmax))/1000))
    set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1],'TickDir','out');
    set(gca,'XTick',1:nmodel,'XTickLabel',models(1:nmodel),'YTick',0:0.2:1);
    xlabel('model');
    ylabel('model probability');
    fig = helper_ILL(gcf,pbar,HGT,false);
    
    set(fig,'PaperPositionMode','manual', ...
        'PaperPosition',[2.5,13,16,HGT],'PaperUnits','centimeters', ...
        'PaperType','A4','PaperOrientation','portrait');  

if isfield(cfg,'saveornot') && isfield(cfg,'figname') && isfield(cfg,'figpath')
    if cfg.saveornot
        print(fullfile(cfg.figpath,cfg.figname),'-painters','-dpdf');
        fprintf('\nfigure %s saved in %s !\n',cfg.figname,cfg.figpath)
    end
end
    
    
end