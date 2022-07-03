function [c] = getc(t,s,r,w,getmat)
%  GETC  Get reversal/repetition curves

% check input arguments
if nargin < 5
    getmat = false;
end
if nargin < 4
    w = [];
end
if nargin < 3
    error('missing input arguments!');
end

ndat = size(r,2); % number of datasets

% find reversals
irev = 1+find(s(2:end) ~= s(1:end-1) & t(2:end) > 1); % indices of reversals
ioff = -4:+3; % index offsets around reversals
nrev = numel(irev); % number of reversals

% create output structure
c = [];
c.ndat = ndat; % number of datasets
c.nrev = nrev; % number of reversals

% compute reversal/repetition curves
crev = nan(8,ndat,nrev); % reversal curve
crep = nan(8,ndat,nrev); % repetition curve
for i = 1:nrev
    crev(:,:,i) = r(irev(i)+ioff,:)*s(irev(i)) == 1;
    crep(:,:,i) = r(irev(i)+ioff,:) == r(irev(i)+ioff-1,:);
end
crev = mean(crev,3);
crep = mean(crep,3);
c.rev_all = crev;
c.rev_avg = mean(crev,2);
c.rev_std = std(crev,[],2);
c.rep_all = crep;
c.rep_avg = mean(crep,2);
c.rep_std = std(crep,[],2);

if ~isempty(w)
    % compute reversal/repetition curves wrt stimulus weight
    if ~all(w > 0)
        error('stimulus weights should be positive!');
    end
    nw = numel(unique(w));
    c.rev_w_all = nan(8,nw,ndat);
    c.rev_w_avg = nan(8,nw);
    c.rev_w_std = nan(8,nw);
    c.rep_w_all = nan(8,nw,ndat);
    c.rep_w_avg = nan(8,nw);
    c.rep_w_std = nan(8,nw);
    for iw = 1:nw
        jrev = intersect(irev,find(w == iw));
        nrev = numel(jrev);
        crev = nan(8,ndat,nrev);
        crep = nan(8,ndat,nrev);
        for i = 1:nrev
            crev(:,:,i) = r(jrev(i)+ioff,:)*s(jrev(i)) == 1;
            crep(:,:,i) = r(jrev(i)+ioff,:) == r(jrev(i)+ioff-1,:);
        end
        crev = mean(crev,3);
        crep = mean(crep,3);
        c.rev_w_all(:,iw,:) = reshape(crev,[8,1,ndat]);
        c.rev_w_avg(:,iw) = mean(crev,2);
        c.rev_w_std(:,iw) = std(crev,[],2);
        c.rep_w_all(:,iw,:) = reshape(crep,[8,1,ndat]);
        c.rep_w_avg(:,iw) = mean(crep,2);
        c.rep_w_std(:,iw) = std(crep,[],2);
    end
end

if getmat
    % compute consistency matrix
    cmat = nan(8,8,ndat,nrev);
    for i = 1:nrev
        for j = 1:8
            for k = j+1:8 % upper triangular part
                cmat(j,k,:,i) = r(irev(i)+ioff(j),:) == r(irev(i)+ioff(k),:);
            end
        end
    end
    cmat = mean(cmat,4);
    c.mat_all = cmat;
    c.mat_avg = mean(cmat,3);
    c.mat_std = std(cmat,[],3);
end

end