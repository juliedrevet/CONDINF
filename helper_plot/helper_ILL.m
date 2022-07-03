function f = helper_ILL(f,pbar,h,mkr_onoff,alpha_onoff)

if nargin < 2
    pbar = 1;
    h = 5;
    mkr_onoff = true;
    mkr_siz = 5;
    alpha_onoff = false;
elseif nargin < 3
    h = 5;
    mkr_onoff = true;
    mkr_siz = 5;
    alpha_onoff = false;
elseif nargin < 4
    mkr_onoff = true;
    alpha_onoff = false;
    if h<5
        mkr_siz = 4.5;
    else
        mkr_siz = 5;
    end
elseif nargin < 5
    mkr_siz = mkr_onoff;
    mkr_onoff = true;
    alpha_onoff = false;
elseif nargin == 5
    mkr_siz = mkr_onoff;
    mkr_onoff = true;
end

set(f,'units','centimeters','innerposition',[0 0 16 h]);

txt_siz = 7.2;
fnt_siz = 8;


% change all axes from figure f to color black
axes = findobj(f, 'type', 'axes');
for a = 1:length(axes)
    if axes(a).YColor <= [1 1 1]
        axes(a).YColor = [0 0 0];
    end
    if axes(a).XColor <= [1 1 1]
        axes(a).XColor = [0 0 0];
    end
end

% access Axes children
axes = findobj(f, 'type', 'axes');
for a = 1:length(axes)
    for ic = 1:numel(axes(a).Children)
        % annotation Fontsize and FontName
        if strcmp(axes(a).Children(ic).Type,'text')
            axes(a).Children(ic).FontSize = txt_siz;
            axes(a).Children(ic).FontName = 'Helvetica';
        end
        % markers 'o' size (from scatter or plot)
        if mkr_onoff
            try strcmp(axes(a).Children(ic).Marker,'o'); % couldn't use isfield, strange...
                try axes(a).Children(ic).MarkerSize = mkr_siz;
                catch;end
                try  axes(a).Children(ic).SizeData= mkr_siz*10;
                catch;end
            catch
            end
        end
        if ~alpha_onoff
            % remove alpha transparency
            try    axes(a).Children(ic).FaceAlpha = 1; % for patch or fill
            catch;end
            try    axes(a).Children(ic).MarkerFaceAlpha = 1; % for scatter
            catch;end
            try    axes(a).Children(ic).MarkerEdgeAlpha = 1; % for scatter
            catch;end
        end
    end
    try axes(a).Legend.FontSize = txt_siz;catch;end
%     % remove title
%     axes(a).Title = [];
%     % remove \newline from XtickLabels
%     axes(a).XTickLabel(contains(axes(a).XTickLabel,'\newline'))=replace(axes(a).XTickLabel(contains(axes(a).XTickLabel,'\newline')),'\newline',' ');
%      % remove \newline from YtickLabels
%      axes(a).YTickLabel(contains(axes(a).YTickLabel,'\newline'))=replace(axes(a).YTickLabel(contains(axes(a).YTickLabel,'\newline')),'\newline',' ');
%      % remove suptitles
%      if strcmp(axes(a).Tag,'suptitle')
%          cla(axes(a));
%      end
end
% change general font, plotboxaspectratio, etc...
axes = findobj(f, 'type', 'axes');
for a = 1:length(axes)
    set(axes(a),'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1]);
    set(axes(a),'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1));
    set(axes(a),'FontName','Helvetica','FontSize',txt_siz);
    axes(a).XLabel.FontSize = fnt_siz;
    axes(a).YLabel.FontSize = fnt_siz;
end

set(f,'units','centimeters','innerposition',[0 0 16 h]);
set(f,'PaperPositionMode','manual', ...
            'PaperPosition',[2.5,13,16,h],'PaperUnits','centimeters', ...
            'PaperType','A4','PaperOrientation','portrait');
end