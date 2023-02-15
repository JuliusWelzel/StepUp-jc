function freqbands = getBounds(data,srate,fbounds,fois)
% getBounds get frequency boundaries from multichannel time series recording
% based in Cohen 2020
%
% INPUT:
%       data -> multichannel times series
%       srate -> sampling rate
%       fbounds -> Frequency boundaries (upper and lower limit)


nbchan  = min(size(data));
dipfrex = fois;
pnts    = length(data);


%% frex params

numfrex = diff(fbounds);
frex = logspace(log10(fbounds(1)),log10(fbounds(2)),numfrex);
stds = linspace(1,1,numfrex);

onsets = srate*2:2*srate:pnts-srate*2;
snipn = 2*srate;

% initialize
[evals,evecs,maps] = deal(zeros(numfrex,nbchan));

%% create R

% full R
R = zeros(nbchan,nbchan);
for segi=1:length(onsets)
    snipdat = data(:,onsets(segi):onsets(segi)+snipn);
    snipdat = bsxfun(@minus,snipdat,mean(snipdat,2));
    R = R + snipdat*snipdat'/snipn;
end
R = R/segi;

% regularized R
gamma = .01;
Rr = R*(1-gamma) + eye(nbchan)*gamma*mean(eig(R));

%% loop over frequencies

for fi=1:numfrex
    
    % filter data
    fdat = filterFGx(data,srate,frex(fi),stds(fi));
    
    %%% compute S
    % full S
    S = zeros(nbchan,nbchan);
    for segi=1:length(onsets)
        snipdat = fdat(:,onsets(segi):onsets(segi)+snipn);
        snipdat = bsxfun(@minus,snipdat,mean(snipdat,2));
        S = S + snipdat*snipdat'/snipn;
    end
    % global variance normalize (optional; this scales the eigenspectrum)
    S = S / (std(S(:))/std(R(:)));

    % GED
    [W,L] = eig(S,Rr);
    [evals(fi,:),sidx] = sort(diag(L),'descend');
    W = W(:,sidx);
    
    % store top component map and eigenvector
    maps(fi,:) = W(:,1)'*S;
    evecs(fi,:) = W(:,1);
    
end

%% correlation matrix for clustering

E = zscore(evecs,[],2);
evecCorMat = (E*E'/(nbchan-1)).^2;

%%

figure(2), clf, set(gcf,'color','w')
contourf(frex,frex,1-evecCorMat,40,'linecolor','none'), hold on
% set(gca,'clim',[0 1],'xscale','log','yscale','log','xtick',round(logspace(log10(1),log10(numfrex),14),1),'ytick',round(logspace(log10(1),log10(numfrex),14),1),'fontsize',15)
set(gca,'clim',[0 1])
xlabel('Frequency (Hz)'), ylabel('Frequency (Hz)') 
axis square, axis xy, colormap bone
title('Eigenvectors correlation matrix')

% box
for i=1:length(dipfrex)
    tbnds = dipfrex{i}; dsearchn(frex',dipfrex{i}');
    plot(tbnds,[1 1]*tbnds(1),'r--','linew',2)
    plot(tbnds,[1 1]*tbnds(2),'r--','linew',2)
    plot([1 1]*tbnds(1),tbnds,'r--','linew',2)
    plot([1 1]*tbnds(2),tbnds,'r--','linew',2)
end

h = colorbar;
set(h,'ticklabels',{'1','.8','.6','.4','.2'},'Ticks',0:.2:1,'fontsize',15)

        
% determine the optimal epsilon value

% range of epsilon parameters
nepsis  = 50;
epsis   = linspace(.001,0.2,nepsis);
qvec    = nan(nepsis,1);

for thi=1:length(epsis)
    
    % scan
    freqbands = dbscan(evecCorMat,epsis(thi),3,'Distance','Correlation');
    
    % compute q
    qtmp = zeros(max(freqbands),1);
    MA = false(size(evecCorMat));
    for i=1:max(freqbands)
        M = false(size(evecCorMat));
        M(freqbands==i,freqbands==i) = 1;
        qtmp(i) = mean(mean(evecCorMat(M))) / mean(mean(evecCorMat(~M)));
        MA = MA+M;
    end
    qvec(thi) = mean(qtmp) + log(mean(MA(:)));
end

% run it again on the best epsilon value
[~,epsiidx] = findpeaks(qvec,'NPeaks',1,'SortStr','descend');
if isempty(epsiidx), epsiidx = round(nepsis/2); end
freqbands = dbscan(evecCorMat,epsis(epsiidx),3,'Distance','Correlation');

% draw empirical bounds on correlation map

for i=1:max(freqbands)
    
    tbnds = frex(freqbands==i);
    tbnds = tbnds([1 end]);
    
    % box
    plot(tbnds,[1 1]*tbnds(1),'m','linew',2)
    plot(tbnds,[1 1]*tbnds(2),'m','linew',2)
    plot([1 1]*tbnds(1),tbnds,'m','linew',2)
    plot([1 1]*tbnds(2),tbnds,'m','linew',2)
end


freqbands = tbnds;
end

