function modelweights = DempsterFor3models(Data)
%DEMPSTERFOR3MODELS.m uses Dempster-Shafer Theory to determine the weight 
%for surrogate model combinations with three models
%--------------------------------------------------------------------------
%Copyright (c) 2013 by Juliane Mueller
%
%This file is part of MATSuMoTo.m - the MATLAB Surrogate Model Toolbox
%MATSuMoTo is free software: you can redistribute it and/or modify it under
%the terms of the GNU General Public License as published by the Free 
%Software Foundation, either version 3 of the License, or (at your option) 
%any later version.
%
%MATSuMoTo is distributed in the hope that it will be useful, but WITHOUT 
%ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
%FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for 
%more details.
%
%You should have received a copy of the GNU General Public License along 
%with MATSuMoTo.  If not, see <http://www.gnu.org/licenses/>.
%--------------------------------------------------------------------------
%Author information
%Juliane Mueller
%juliane.mueller2901@gmail.com
%--------------------------------------------------------------------------
%
%Input:
%Data - structure with all problem information
%Output:
%modelweights - vector with weights for each model
%--------------------------------------------------------------------------

m=size(Data.S,1); %numer of already sampled points
yPred_all=zeros(m,3); %initialize matrix for predicted function values

%% determine whether k-fold or leave-one-out cross-validation in needed
if m>50  %k-fold cross validation
    multiplier=ceil(m/50)-1;
    kfold=10*multiplier;
    forcekfold=true;
    mix=randperm(m);
else %leave-one-out cross-validation
    forcekfold=false;
end 

Data.numberOfModels=3; %number of models in ensemble
if strcmp(Data.surrogate_model,'MIX_RcPcM') %ensemble of cubic RBF, full cubic polynomial, and MARS
    if ~forcekfold  %leave-one-out cross-validation
        for jj=1:m
            %fit cubic RBF model to training set
            [lambda, gamma]=RBF([Data.S(1:jj-1,:);Data.S(jj+1:end,:)],...
                [Data.Ymed(1:jj-1,:);Data.Ymed(jj+1:end,:)],'cubic'); 
            %re-predict function value for validation set with cubic RBF
            yPred_all(jj,1)=RBF_eval(Data.S(jj,:),[Data.S(1:jj-1,:);...
                Data.S(jj+1:end,:)],lambda,gamma,'cubic');
            % fit full cubic polynomial to training set
            beta=POLY([Data.S(1:jj-1,:);Data.S(jj+1:end,:)],...
                [Data.Ymed(1:jj-1,:);Data.Ymed(jj+1:end,:)],'cub'); 
            %re-predict function value for validation set with full cubic polynomial
            yPred_all(jj,2)=POLY_eval(Data.S(jj,:),beta,'cub');  
            % fit MARS model to training set
            mmodel = aresbuild([Data.S(1:jj-1,:);Data.S(jj+1:end,:)],...
                [Data.Ymed(1:jj-1,:);Data.Ymed(jj+1:end,:)]); 
            %re-predict function value for validation set with MARS
            yPred_all(jj,3)=arespredict(mmodel,Data.S(jj,:)); 
        end      
    else
        Ss=Data.S(mix,:); %resorted sample sites
        Ys=Data.Ymed(mix);
        for jj=1:ceil(m/kfold)
            %validation set
            if jj*kfold<=m
                validation_S=Ss((jj-1)*kfold+1:jj*kfold,:);
            else %left-overs
                validation_S=Ss((jj-1)*kfold+1:end,:);
            end
            % training set
            if jj ==1 %first group
                trainset_S=Ss(jj*kfold+1:end,:);
                trainset_Y=Ys(jj*kfold+1:end,:);
            elseif 2<=jj && jj<= floor(m/kfold)  && mod(m,kfold) >0 ||...
                    2<=jj && mod(m,kfold)==0 && jj*kfold<m
                    %group 2 till last full group %mod(m,kfold)~=0
                trainset_S=[Ss(1:(jj-1)*kfold,:);Ss(jj*kfold+1:end,:)];
                trainset_Y=[Ys(1:(jj-1)*kfold,:);Ys(jj*kfold+1:end,:)];
            else %left-overs
                trainset_S=Ss(1:(jj-1)*kfold,:);
                trainset_Y=Ys(1:(jj-1)*kfold,:);
            end
            if jj ==1 || 1<=jj && jj<= floor(m/kfold)  && mod(m,kfold) >0 ||...
                    2<=jj && mod(m,kfold)==0 && jj*kfold<m
                %fit cubic RBF model to training set
                [lambda, gamma]=RBF(trainset_S,trainset_Y,'cubic');
                %re-predict function value for validation set with cubic RBF
                yPred_all((jj-1)*kfold+1:jj*kfold,1)=RBF_eval(validation_S,...
                    trainset_S,lambda,gamma,'cubic');
                % fit full cubic polynomial to training set
                beta=POLY(trainset_S,trainset_Y,'cub'); 
                %re-predict function value for validation set with full cubic polynomial
                yPred_all((jj-1)*kfold+1:jj*kfold,2)=POLY_eval(validation_S,beta,'cub');
                % fit MARS model to training set
                mmodel = aresbuild(trainset_S,trainset_Y); 
                %re-predict function value for validation set with MARS
                yPred_all((jj-1)*kfold+1:jj*kfold,3)=arespredict(mmodel,validation_S);
            else %left-overs
                id_a=(jj-1)*kfold+1;
                if mod(m,kfold)>0
                    id_e=(jj-1)*kfold+mod(m,kfold);
                else
                    id_e=jj*kfold;
                end
                %fit cubic RBF model to training set
                [lambda, gamma]=RBF(trainset_S,trainset_Y,'cubic');
                %re-predict function value for validation set with cubic RBF
                yPred_all(id_a:id_e,1)=RBF_eval(validation_S,trainset_S,...
                    lambda,gamma,'cubic');
                % fit full cubic polynomial to training set
                beta=POLY(trainset_S,trainset_Y,'cub'); 
                %re-predict function value for validation set with full cubic polynomial
                yPred_all(id_a:id_e,2)=POLY_eval(validation_S,beta,'cub');
                % fit MARS model to training set
                mmodel = aresbuild(trainset_S,trainset_Y); 
                %re-predict function value for validation set with MARS
                yPred_all(id_a:id_e,3)=arespredict(mmodel,validation_S);
                yPred_all=sortrows([yPred_all,mix(:)],Data.numberOfModels+1);
                yPred_all=yPred_all(:,1:Data.numberOfModels);
            end
        end  
    end
%add further 3-model combinations here     
%elseif strcmp(Data.surrogate_model,'YourModel')
end
        
%% characteristics of individual models 
Corcoef=cc_calc(Data.Ymed,yPred_all); %correlation coefficient
RMSE=RMSE_calc(Data.Ymed,yPred_all); %root mean squared error
MAE=MAE_calc(Data.Ymed,yPred_all); %maximum absolute errors
MAD=MAD_calc(Data.Ymed,yPred_all); %median absolute deviation
    
%correlation coefficients    
Corcoef_t=Corcoef;
if all(Corcoef_t==Corcoef_t(2)) %same CC value for wach model
    Corcoef_t=1/Data.numberOfModels*ones(Data.numberOfModels,1);
end
if any(Corcoef_t<0) %scale correlation coefficients to interval 0,1
    Corcoef_t=(Corcoef_t-min(Corcoef_t))./max(Corcoef_t-min(Corcoef_t));
end
if all(isnan(Corcoef))
    Corcoef_t=1/Data.numberOfModels*ones(Data.numberOfModels,1);
else
    Corcoef_t(Corcoef_t==0)=0.000001;
    Corcoef_t=Corcoef_t(:)./(sum(Corcoef_t(:)));
end
%Root mean squared errors
if all(RMSE==0)
    RMSE_t=1/Data.numberOfModels*ones(Data.numberOfModels,1);
else
    RMSE_t=(1./RMSE(:))./sum(1./RMSE(:));
end
% Maximum absolute errors
if all(MAE==0)
    MAE_t=1/Data.numberOfModels*ones(Data.numberOfModels,1);
else
    MAE_t=(1./MAE(:))./sum(1./MAE(:));
end
%Median absolute deviation
if all(MAD==0)
    MAD_t=1/Data.numberOfModels*ones(Data.numberOfModels,1);
else
    MAD_t=(1./MAD(:))./(sum(1./MAD));
end

%combine model information
comb=[1 NaN NaN; NaN 2 NaN; NaN NaN 3]; %individual models
OVpcr=[comb Corcoef_t(:), RMSE_t(:), MAE_t(:), MAD_t(:)]; %evidences of individual models
    
intround=true;
%DST to combine evidences
for yy =2:Data.numberOfModels
    if intround
        compD=OVpcr(:,1:Data.numberOfModels+1);
        intround=false;
    else
        compD=Intersections;
    end
    OVnew=[comb,OVpcr(:,Data.numberOfModels+yy)];%Transform data into evidences
    Intersections=dempster_rule(compD,OVnew,Data);%combination rule
end
%calculate believes, plausibilities and pignistic probabilities
Intervals=Dempster_belpl(Intersections,Data);
%predict function values with all possible 2-model combinations and with 3-model combination    
Ymixtest=model2combi2_3mod(Intervals,yPred_all); 
%characteristics of mixture
Corcoef=[Corcoef;cc_calc(Data.Ymed,Ymixtest)];
if all(Corcoef_t==Corcoef_t(2)) %same CC value for wach model
    Corcoef=1/length(Corcoef)*ones(length(Corcoef),1);
end
if any(Corcoef<0)
    Corcoef=(Corcoef-min(Corcoef))./max(Corcoef-min(Corcoef));
end
Corcoef(Corcoef==0)=0.000001;
RMSE=[RMSE;RMSE_calc(Data.Ymed,Ymixtest)]; %root mean squared errors
if all(RMSE==0)  %Same RMSE for each model
    RMSE=1/length(RMSE)*ones(length(RMSE),1);
end
MAE=[MAE;MAE_calc(Data.Ymed,Ymixtest)]; %maximum absolute errors
if all(MAE==0)  %Same MAE for each model
    MAE=1/length(MAE)*ones(length(MAE),1);
end
MAD=[MAD;MAD_calc(Data.Ymed,Ymixtest)]; % median absolute deviation
if all(MAD==0)  %Same MAD for each model
    MAD=1/length(MAD)*ones(length(MAD),1);
end

%% Combine pure models to mixed model     
%all available model combinations
comb=[1 NaN NaN; NaN 2 NaN; NaN NaN 3; 1 2 NaN; 1 NaN 3; NaN 2 3; 1 2 3]; 

%matrix containing model combination and corresponding evidences (high numbers mean high evidence)
OVpcr=[comb Corcoef(:)./(sum(Corcoef(:))) (1./RMSE(:))./sum(1./RMSE(:))...
    (1./MAE(:))./sum(1./MAE(:)) (1./MAD(:))./(sum(1./MAD))];
clear Corcoef; clear RMSE; clear MAE; clear MAD
clear Corcoef_t; clear RMSE_t; clear MAE_t; clear MAD_t;

intround=true;
%DST
for yy =2:Data.numberOfModels
    if intround
        compD=OVpcr(:,1:Data.numberOfModels+1);
        intround=false;
    else
        compD=Intersections;
    end
    OVnew=[comb,OVpcr(:,Data.numberOfModels+yy)]; %Transform data into evidences
    Intersections=dempster_rule(compD,OVnew,Data); %combination rule
end

%calculate believes, plausibilities and pignistic probabilities
Intervals=Dempster_belpl(Intersections,Data);
newmodel=Intervals(1,1:Data.numberOfModels);
newmodel(isnan(newmodel))=[];
%determine the weights of the models in the combination
w=weights_in_combi(Intervals,newmodel,Data);
modelweights=zeros(1,Data.numberOfModels);
modelweights(newmodel)=w/sum(w);
end%function