function modelweights = DempsterFor2models(Data)
%DEMPSTERFOR2MODELS.m uses Dempster-Shafer Theory to determine the weight 
%for surrogate model combinations with two models
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
%input: 
%Data - structure  with all model information obtained so far
%Output:
%modelweights - vector with weight of each model
%--------------------------------------------------------------------------

m=size(Data.S,1); %number of points already sampled
yPred_all=zeros(m,2);  %initialize array for objective function value predictions

%% determine whether k-fold or leave-one-out cross-validation
if m>50  %k-fold cross validation
    multiplier=ceil(m/50)-1;
    kfold=10*multiplier;
    forcekfold=true;
    mix=randperm(m);
else %leave-one-out cross-validation
    forcekfold=false;
end    

Data.numberOfModels=2; %number of models in the ensemble

%% cross-validation for mixtures
%Ensemble: cubic RBF & MARS    
if strcmp(Data.surrogate_model,'MIX_RcM')  
    if ~forcekfold  % leave-1-out cross-validation
        for jj=1:m
            %fit cubic RBF model to training set
            [lambda, gamma]=RBF([Data.S(1:jj-1,:);Data.S(jj+1:end,:)],...
                [Data.Ymed(1:jj-1,:);Data.Ymed(jj+1:end,:)],'cubic');
            %re-predict function value for validation set with cubic RBF
            yPred_all(jj,1)=RBF_eval(Data.S(jj,:),[Data.S(1:jj-1,:);...
                Data.S(jj+1:end,:)],lambda,gamma,'cubic');
            % fit MARS model to training set
            mmodel = aresbuild([Data.S(1:jj-1,:);Data.S(jj+1:end,:)],...
                [Data.Ymed(1:jj-1,:);Data.Ymed(jj+1:end,:)]); 
            %re-predict function value for validation set with MARS
            yPred_all(jj,2)=arespredict(mmodel,Data.S(jj,:)); 
        end
    else %k-fold cross-validation
        Ss=Data.S(mix,:); %randomly sorted sample sites
        Ys=Data.Ymed(mix); %randomly sorted function values
        for jj=1:ceil(m/kfold)
            %determine validation set
            if jj*kfold<=m
                validation_S=Ss((jj-1)*kfold+1:jj*kfold,:); %validation points
            else %left-overs
                validation_S=Ss((jj-1)*kfold+1:end,:); %validation points
            end
            %determine training set
            if jj ==1 %first group
                trainset_S=Ss(jj*kfold+1:end,:); %training points
                trainset_Y=Ys(jj*kfold+1:end,:); %function values of training points
            elseif 2<=jj && jj<= floor(m/kfold)  && mod(m,kfold) >0 ||...
                    2<=jj && mod(m,kfold)==0 && jj*kfold<m
                    %group 2 till last full group %mod(m,kfold)~=0
                trainset_S=[Ss(1:(jj-1)*kfold,:);Ss(jj*kfold+1:end,:)]; %training points
                trainset_Y=[Ys(1:(jj-1)*kfold,:);Ys(jj*kfold+1:end,:)]; %function values of training points
            else %left-overs
                trainset_S=Ss(1:(jj-1)*kfold,:); %training points
                trainset_Y=Ys(1:(jj-1)*kfold,:); %function values of training points
            end
            if jj==1 || 1<=jj && jj<= floor(m/kfold)  && mod(m,kfold) >0 ||...
                    2<=jj && mod(m,kfold)==0 && jj*kfold<m
                %fit cubic RBF model to training set
                [lambda, gamma]=RBF(trainset_S,trainset_Y,'cubic');
                %re-predict function value for validation set with cubic RBF
                yPred_all((jj-1)*kfold+1:jj*kfold,1)=RBF_eval(validation_S,...
                    trainset_S,lambda,gamma,'cubic');
                % fit MARS model to training set
                mmodel = aresbuild(trainset_S,trainset_Y);
                %re-predict function value for validation set with MARS
                yPred_all((jj-1)*kfold+1:jj*kfold,2)=arespredict(mmodel,validation_S); 
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
                % fit MARS model to training set
                mmodel = aresbuild(trainset_S,trainset_Y); 
                %re-predict function value for validation set with MARS
                yPred_all(id_a:id_e,2)=arespredict(mmodel,validation_S); 
                yPred_all=sortrows([yPred_all,mix(:)],Data.numberOfModels+1);
                yPred_all=yPred_all(:,1:Data.numberOfModels);
            end
        end 
    end
%Ensemble: cubic RBF & full cubic polynomial    
elseif strcmp(Data.surrogate_model,'MIX_RcPc')  
    if ~forcekfold  % leave-one-out cross-validation
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
        end
    else %k-fold cross-validation
        Ss=Data.S(mix,:); %resorted sample sites
        Ys=Data.Ymed(mix);
        for jj=1:ceil(m/kfold)
            %validation set
            if jj*kfold<=m
                validation_S=Ss((jj-1)*kfold+1:jj*kfold,:);
            else %left-overs
                validation_S=Ss((jj-1)*kfold+1:end,:);
            end
            %training set
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
                %fit full cubic polynomial to training set
                beta=POLY(trainset_S,trainset_Y,'cub'); 
                %re-predict function value for validation set with full cubic polynomial
                yPred_all(id_a:id_e,2)=POLY_eval(validation_S,beta,'cub');
                yPred_all=sortrows([yPred_all,mix(:)],Data.numberOfModels+1);
                yPred_all=yPred_all(:,1:Data.numberOfModels);
            end
        end 
    end
%Ensemble: cubic RBF & full quadratic polynomial    
elseif strcmp(Data.surrogate_model,'MIX_RcPq')  
    if ~forcekfold  % leave-one-out cross-validation
        for jj=1:m
            %fit cubic RBF model to training set
            [lambda, gamma]=RBF([Data.S(1:jj-1,:);Data.S(jj+1:end,:)],...
                [Data.Ymed(1:jj-1,:);Data.Ymed(jj+1:end,:)],'cubic'); 
            %re-predict function value for validation set with cubic RBF
            yPred_all(jj,1)=RBF_eval(Data.S(jj,:),[Data.S(1:jj-1,:);...
                Data.S(jj+1:end,:)],lambda,gamma,'cubic');
            % fit full quadratic polynomial to training set
            beta=POLY([Data.S(1:jj-1,:);Data.S(jj+1:end,:)],...
                [Data.Ymed(1:jj-1,:);Data.Ymed(jj+1:end,:)],'quad'); 
            %re-predict function value for validation set with full quadratic polynomial
            yPred_all(jj,2)=POLY_eval(Data.S(jj,:),beta,'quad');
        end
    else %k-fold cross-validation
        Ss=Data.S(mix,:); %resorted sample sites
        Ys=Data.Ymed(mix);
        for jj=1:ceil(m/kfold)
            %validation set
            if jj*kfold<=m
                validation_S=Ss((jj-1)*kfold+1:jj*kfold,:);
            else %left-overs
                validation_S=Ss((jj-1)*kfold+1:end,:);
            end
            %training set
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
            if jj ==1 ||1<=jj && jj<= floor(m/kfold)  && mod(m,kfold) >0 ||...
                    2<=jj && mod(m,kfold)==0 && jj*kfold<m
                %fit cubic RBF model to training set
                [lambda, gamma]=RBF(trainset_S,trainset_Y,'cubic');
                %re-predict function value for validation set with cubic RBF
                yPred_all((jj-1)*kfold+1:jj*kfold,1)=RBF_eval(validation_S,...
                    trainset_S,lambda,gamma,'cubic');
                % fit full quadratic polynomial to training set
                beta=POLY(trainset_S,trainset_Y,'quad'); 
                %re-predict function value for validation set with full quadratic polynomial
                yPred_all((jj-1)*kfold+1:jj*kfold,2)=POLY_eval(validation_S,beta,'quad');
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
                % fit full quadratic polynomial to training set
                beta=POLY(trainset_S,trainset_Y,'quad'); 
                %re-predict function value for validation set with full quadratic polynomial
                yPred_all(id_a:id_e,2)=POLY_eval(validation_S,beta,'quad');
                yPred_all=sortrows([yPred_all,mix(:)],Data.numberOfModels+1);
                yPred_all=yPred_all(:,1:Data.numberOfModels);
            end
        end 
    end
%Ensemble: cubic RBF & reduced cubic polynomial
elseif strcmp(Data.surrogate_model,'MIX_RcPcr')
    if ~forcekfold  %leave-one-out cross-validation
        for jj=1:m
            %fit cubic RBF model to training set
            [lambda, gamma]=RBF([Data.S(1:jj-1,:);Data.S(jj+1:end,:)],...
                [Data.Ymed(1:jj-1,:);Data.Ymed(jj+1:end,:)],'cubic'); 
            %re-predict function value for validation set with cubic RBF
            yPred_all(jj,1)=RBF_eval(Data.S(jj,:),[Data.S(1:jj-1,:);...
                Data.S(jj+1:end,:)],lambda,gamma,'cubic'); 
            % fit reduced cubic polynomial to training set
            beta=POLY([Data.S(1:jj-1,:);Data.S(jj+1:end,:)],...
                [Data.Ymed(1:jj-1,:);Data.Ymed(jj+1:end,:)],'cubr'); 
            %re-predict function value for validation set with reduced cubic polynomial
            yPred_all(jj,2)=POLY_eval(Data.S(jj,:),beta,'cubr');
        end
    else %k-fold cross-validation
        Ss=Data.S(mix,:); %resorted sample sites
        Ys=Data.Ymed(mix);
        for jj=1:ceil(m/kfold)
            %validation set
            if jj*kfold <= m
                validation_S=Ss((jj-1)*kfold+1:jj*kfold,:);
            else %left-overs
                validation_S=Ss((jj-1)*kfold+1:end,:);
            end
            %validation set
            if jj ==1 %first group
                trainset_S=Ss(jj*kfold+1:end,:);
                trainset_Y=Ys(jj*kfold+1:end,:);
            elseif  2<=jj && jj<= floor(m/kfold)  && mod(m,kfold) >0 ||...
                2<=jj && mod(m,kfold)==0 && jj*kfold<m
                %group 2 till last full group
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
                % fit reduced cubic polynomial to training set
                beta=POLY(trainset_S,trainset_Y,'cubr'); 
                %re-predict function value for validation set with reduced cubic polynomial
                yPred_all((jj-1)*kfold+1:jj*kfold,2)=POLY_eval(validation_S,beta,'cubr');
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
                % fit reduced cubic polynomial to training set
                beta=POLY(trainset_S,trainset_Y,'cubr'); 
                %re-predict function value for validation set with reduced cubic polynomial
                yPred_all(id_a:id_e,2)=POLY_eval(validation_S,beta,'cubr');
                yPred_all=sortrows([yPred_all,mix(:)],Data.numberOfModels+1);
                yPred_all=yPred_all(:,1:Data.numberOfModels);
            end
        end
    end
%Ensemble: cubic RBF & reduced quadratic polynomial
elseif strcmp(Data.surrogate_model,'MIX_RcPqr')
    if ~forcekfold  %leave-one-out cross-validation
        for jj=1:m
            %fit cubic RBF model to training set
            [lambda, gamma]=RBF([Data.S(1:jj-1,:);Data.S(jj+1:end,:)],...
                [Data.Ymed(1:jj-1,:);Data.Ymed(jj+1:end,:)],'cubic'); 
            %re-predict function value for validation set with cubic RBF
            yPred_all(jj,1)=RBF_eval(Data.S(jj,:),[Data.S(1:jj-1,:);...
                Data.S(jj+1:end,:)],lambda,gamma,'cubic'); 
            % fit reduced quadratic polynomial to training set
            beta=POLY([Data.S(1:jj-1,:);Data.S(jj+1:end,:)],...
                [Data.Ymed(1:jj-1,:);Data.Ymed(jj+1:end,:)],'quadr'); 
            %re-predict function value for validation set with reduced quadratic polynomial
            yPred_all(jj,2)=POLY_eval(Data.S(jj,:),beta,'quadr');
        end
    else %k-fold cross-validation
        Ss=Data.S(mix,:); %resorted sample sites
        Ys=Data.Ymed(mix);
        for jj=1:ceil(m/kfold)
            %validation set
            if jj*kfold <= m
                validation_S=Ss((jj-1)*kfold+1:jj*kfold,:);
            else %left-overs
                validation_S=Ss((jj-1)*kfold+1:end,:);
            end
            %validation set
            if jj ==1 %first group
                trainset_S=Ss(jj*kfold+1:end,:);
                trainset_Y=Ys(jj*kfold+1:end,:);
            elseif 2<=jj && jj<= floor(m/kfold)  && mod(m,kfold) >0 ||...
                2<=jj && mod(m,kfold)==0 && jj*kfold<m
                %group 2 till last full group
                trainset_S=[Ss(1:(jj-1)*kfold,:);Ss(jj*kfold+1:end,:)];
                trainset_Y=[Ys(1:(jj-1)*kfold,:);Ys(jj*kfold+1:end,:)];
            else %left-overs
                trainset_S=Ss(1:(jj-1)*kfold,:);
                trainset_Y=Ys(1:(jj-1)*kfold,:);
            end
            if jj ==1 ||1<=jj && jj<= floor(m/kfold)  && mod(m,kfold) >0 ||...
                2<=jj && mod(m,kfold)==0 && jj*kfold<m
                %fit cubic RBF model to training set
                [lambda, gamma]=RBF(trainset_S,trainset_Y,'cubic');
                %re-predict function value for validation set with cubic RBF
                yPred_all((jj-1)*kfold+1:jj*kfold,1)=RBF_eval(validation_S,...
                    trainset_S,lambda,gamma,'cubic');
                % fit reduced quadratic polynomial to training set
                beta=POLY(trainset_S,trainset_Y,'quadr'); 
                yPred_all((jj-1)*kfold+1:jj*kfold,2)=POLY_eval(validation_S,beta,'quadr');
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
                % fit reduced quadratic polynomial to training set
                beta=POLY(trainset_S,trainset_Y,'quadr'); 
                %re-predict function value for validation set with reduced quadratic polynomial
                yPred_all(id_a:id_e,2)=POLY_eval(validation_S,beta,'quadr');
                yPred_all=sortrows([yPred_all,mix(:)],Data.numberOfModels+1);
                yPred_all=yPred_all(:,1:Data.numberOfModels);
            end
        end
    end    
    
%add other 2-model mixtures here  

end
        
%% characteristics of single models 
Corcoef=cc_calc(Data.Ymed,yPred_all); %correlation coefficients
RMSE=RMSE_calc(Data.Ymed,yPred_all); %root mean squared errors
MAE=MAE_calc(Data.Ymed,yPred_all); %maximum absolute error
MAD=MAD_calc(Data.Ymed,yPred_all); %median absolute deviation
   
%Scale correlation coefficients
Corcoef_t=Corcoef;
if Corcoef_t(1)==Corcoef_t(2) %same CC value for wach model
    Corcoef_t=1/Data.numberOfModels*ones(Data.numberOfModels,1);
end
if any(Corcoef_t<0) %scale correlation coefficients to interval [0,1]
    Corcoef_t=(Corcoef_t-min(Corcoef_t))./max(Corcoef_t-min(Corcoef_t));
end
if all(isnan(Corcoef)) %all models get same CC
    Corcoef_t=1/Data.numberOfModels*ones(Data.numberOfModels,1);
else
    Corcoef_t(Corcoef_t==0)=0.000001; %every CC most be positive
    Corcoef_t=Corcoef_t(:)./(sum(Corcoef_t(:)));
end
%scale root mean squared errors
if all(RMSE==0)  %Same RMSE for each model
    RMSE_t=1/Data.numberOfModels*ones(Data.numberOfModels,1);
else
    RMSE_t=(1./RMSE(:))./sum(1./RMSE(:));
end
%scale Maximum absolute errors
if all(MAE==0) %Same MAE for each model
    MAE_t=1/Data.numberOfModels*ones(Data.numberOfModels,1);
else
    MAE_t=(1./MAE(:))./sum(1./MAE(:));
end
%scale median absolute deviation
if all(MAD==0) %Same MAD for each model
    MAD_t=1/Data.numberOfModels*ones(Data.numberOfModels,1);
else
    MAD_t=(1./MAD(:))./(sum(1./MAD));
end

%compute basic probability assignments (BPA)    
comb=[1 NaN; NaN 2]; %individual models
OVpcr=[comb Corcoef_t(:), RMSE_t(:), MAE_t(:), MAD_t(:)]; %evidences
    
intround=true; 
%DST
for yy =2:Data.numberOfModels
    if intround %first model
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
%compute prediction of mixture model
Ymixtest=zeros(size(yPred_all,1),nchoosek(Data.numberOfModels,2));
for ii = 1:nchoosek(Data.numberOfModels,2)
    Ymixtest(:,ii)=Intervals(1,Data.numberOfModels+1)*yPred_all(:,1)+...
        Intervals(1,Data.numberOfModels+1)*yPred_all(:,2); %mixture model predictions at sample points 
end
%compute model characteristics for mixture model
%compute correlation coefficient of prediction by mixture    
Corcoef=[Corcoef;cc_calc(Data.Ymed,Ymixtest)];
if all(Corcoef == Corcoef(3))
    Corcoef = ones(3,1)/3; %same CC for each model
end
if any(Corcoef<0)
    Corcoef=(Corcoef-min(Corcoef))./max(Corcoef-min(Corcoef));
end
Corcoef(Corcoef==0)=0.000001;

RMSE=[RMSE;RMSE_calc(Data.Ymed,Ymixtest)]; %compute RMSE of mixture
if all(RMSE==0)  %Same RMSE for each model
    RMSE=1/3*ones(3,1);
end
MAE=[MAE;MAE_calc(Data.Ymed,Ymixtest)]; %compute MAE of mixture
if all(MAE==0)  %Same MAE for each model
    MAE=1/3*ones(3,1);
end
MAD=[MAD;MAD_calc(Data.Ymed,Ymixtest)]; % compute MAD of mixture     
if all(MAD==0)  %Same MAD for each model
    MAD=1/3*ones(3,1);
end

%% Combine pure models to mixed model     
%all available model combinations (pure models and 2-model mixture)
comb=[1 NaN; NaN 2; 1 2]; 

%matrix containing model combination and corresponding evidences (high numbers mean high evidence)
OVpcr=[comb Corcoef(:)./(sum(Corcoef(:))) (1./RMSE(:))./sum(1./RMSE(:))...
    (1./MAE(:))./sum(1./MAE(:)) (1./MAD(:))./(sum(1./MAD))];
clear Corcoef; clear RMSE; clear MAE; clear MAD
clear Corcoef_t; clear RMSE_t; clear MAE_t; clear MAD_t;

intround=true;
%DST to determine weights in models
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
w=Intervals([2:3],end);%weights_combi(Intervals,ModelParam.newmodel);
modelweights=w/sum(w); %scale model weights so that w1+w2=1

end % function