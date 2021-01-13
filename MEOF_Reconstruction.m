%------------------------%
%% (1) LOAD DATA
%------------------------%

load('Virtual_Stations.mat')
load('TG_Data.mat')
load('Steric_Height.mat')

%----------------------------------------%
%% (2) Spatial Correlation Patterns % MEOF
%----------------------------------------%

inGoB(1).inGoB = inpolygon(LST(:,1),LST(:,2),[-60 -7 -7 10 10 -20 -60],[0 0 40 50 70 70 0]);
inGoB(2).inGoB = inpolygon(LST(:,1),LST(:,2),[-90 -20 -20 -90 -90],[15 15 60 60 15]);
inGoB(3).inGoB = inpolygon(LST(:,1),LST(:,2),[-90 -20 -20 -90 -90],[15 15 60 60 15]);
inGoB(4).inGoB = inpolygon(LST(:,1),LST(:,2),[-180 -60 -60 -83 -102 -102 -180 -180],[0 0 8.5 8.5 24 60 60 0]);
LST(LST(:,1)<0,1) = LST(LST(:,1)<0,1)+360;
inGoB(5).inGoB = inpolygon(LST(:,1),LST(:,2),[100 260 260 100 100],[-10 -10 60 60 -10]);
LST(LST(:,1)>180,1) = LST(LST(:,1)>180,1)-360;
inGoB(6).inGoB = inpolygon(LST(:,1),LST(:,2),[100 180 180 142 142 100 100],[0 0 36 36 50 50 0]);
inGoB(7).inGoB = inpolygon(LST(:,1),LST(:,2),[60 180 180 60 60],[-60 -60 15 15 -60]);
LST(LST(:,1)<0,1) = LST(LST(:,1)<0,1)+360;
inGoB(8).inGoB = inpolygon(LST(:,1),LST(:,2),[80 300 300 80 80],[-70 -70 0 -20 -70]);
inGoB(9).inGoB = inpolygon(LST(:,1),LST(:,2),[120 300 300 120 120],[0 0 50 50 0]);
LST(LST(:,1)>180,1) = LST(LST(:,1)>180,1)-360;

s = find(t>=1980);
for j = 1:9;
    STGoB = detrend(permute(ST(inGoB(j).inGoB,:,:),[2 1 3]),'constant');
    LSTGoB = LST(inGoB(j).inGoB,:);
    ZGoB = Z(inGoB(j).inGoB);
    KST = NaN(size(STGoB,2),1);
    KSTSig = NaN(size(STGoB,2),1);
    KSTP = NaN(size(STGoB,2),1);
    % Spatial Correlation Analysis
    for i = 1:size(STGoB,2);
        isn = find(~isnan(STGoB(:,i)));
        if length(isn)==size(STGoB,1);
            [r,rP] = corrcoef(detrend(STGoB(s,i)),detrend(SIRes(s,j)));
            KST(i,1) = r(2);
            KSTP(i,1) = rP(2);
        else
            KST(i,1) = NaN;
            KSTP(i,1) = NaN;
        end
    end
    s2 = find(KST>0&KSTP<0.05&ZGoB>=300);
    for i = 1:size(MResens,3);
        % MEOF Reconstruction
        smEOF = randi([3 6],1,1); % random smoothing
        texp = randi([92 98],1,1); % random number of EOFs included in the regression
        Pred = movmean(detrend([STGoB(:,s2,1),MResens(:,index==j,i)]),smEOF);      
        [COEFF,SCORE,~,~,explained,~] = pca(Pred);  
        seof = find(cumsum(explained)<=texp);
        NumEOF(i,j) = size(seof,1);
        m = STGoB(:,s2,1)/COEFF(1:size(s2,1),seof)';
        RecTG(j).RECTG(:,:,i) = m(:,seof)*COEFF(size(s2,1)+1:end,seof)';
        EOFRec(:,i,j) = mean(RecTG(j).RECTG(:,:,i),2)+SIATM(:,j);
    end
    SText(j).L = LSTGoB(s2,:);
    SText(j).K = KST(s2,:);SText(j).ST = STGoB(:,s2,:);
end