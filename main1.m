close all; clc;clear
addpath('funs/')
addpath('results/')
addpath('data/')
load('AR_28x20.mat')
X   = X';
gnd = X_label';
X   = double(X);
fea = NormalizeFea(X',0);
% options = [];
% options.ReducedDim = 800;
% [eigvector,eigvalue,meanData,new_data] = PCA(fea',options);
k=68;
for lam1 =  [1e-5,2*1e-5, 5*1e-5,8*1e-5,1e-4,2*1e-4, 5*1e-4,8*1e-4,1e-3,2*1e-3,5*1e-3,8*1e-3,1e-2,2*1e-2,5*1e-2,8*1e-2,1e-1,0.2,0.5,0.8,1,2,5,8,10,50,100]
    for lam2 =0.1% [1e-5,2*1e-5, 5*1e-5,8*1e-5,1e-4,2*1e-4, 5*1e-4,8*1e-4,1e-3,2*1e-3,5*1e-3,8*1e-3,1e-2,2*1e-2,5*1e-2,8*1e-2,1e-1,0.2,0.5,0.8,1,2,5,8,10,50,100]
r=1;
c=max(gnd);
[Z,F] = LapdoubleL21SR(fea,lam1,lam2,c-1);
% [Z] = L21SR(fea,lambda,r);
% [Z1] = doubleL21SR(fea,lam1);
% [Z2] = LRR(new_data',lambda);
% [Z,E] = lrra(X,X,lambda,1);
% figure;imshow(Z,[])
% figure;show((fea*Z)')
% figure;plot(Z(:,11),'b')


% L1 = BuildAdjacency(thrC(Z1),1);
L = BuildAdjacency(thrC(Z),1);
for i = 1:10
% P_label1{i} = SC(L1,c);
% label1{i} = reshape(P_label1{i},1,[]);
% res1{i} = bestMap(gnd,label1{i});
% acc1(i) = length(find(gnd == res1{i}))/length(gnd);
% MIhat1(i) = MutualInfo(gnd,res1{i});


P_label{i} = SC(L,c);
label{i} = reshape(P_label{i},1,[]);
res{i} = bestMap(gnd,label{i});
confmatrix{i} = cfmatrix(gnd, res{i}, 1:c, 0);
TP{i}      = diag(confmatrix{i});
FP{i}      = sum(confmatrix{i});
FP{i}      = FP{i}(:)-TP{i}(:);
FN{i}      = sum(confmatrix{i},2);
FN{i}      = FN{i}(:)-TP{i}(:);
TN{i}      = sum(sum(confmatrix{i}))-(TP{i}+FP{i}+FN{i});
RI(i) = sum((TP{i}+TN{i})./(TP{i}+FP{i}+TN{i}+FN{i}))/c;
acc(i) = sum(TP{i}./(TP{i}+FP{i}+TN{i}+FN{i}));
Precision{i} = TP{i}./(TP{i}+FP{i});
Recall{i} = TP{i}./(TP{i}+FN{i});
F1m     = (2*Precision{i}.*Recall{i})./(Precision{i}+Recall{i});
F1m(isnan(F1m)) = 0;
F1(i)   = sum(F1m)/c;
MIhat(i) = MutualInfo(gnd,res{i});
TPs(i)      = sum(TP{i})/c;
FPs(i)      = sum(FP{i})/c;
TNs(i)      = sum(TN{i})/c;
FNs(i)      = sum(FN{i})/c;
Precisions(i)= sum(Precision{i})/c;
Recalls(i)   = sum(Recall{i})/c;
%% F predict
MAXiter = 500;
REPlic = 20;
P_labelF{i} =  litekmeans(F,c,'MaxIter',MAXiter,'Replicates',REPlic);
labelF{i} = reshape(P_labelF{i},1,[]);
resF{i} = bestMap(gnd,labelF{i});
confmatrixF{i} = cfmatrix(gnd, resF{i}, 1:c, 0);
TPF{i}      = diag(confmatrixF{i});
FPF{i}      = sum(confmatrixF{i});
FPF{i}      = FPF{i}(:)-TPF{i}(:);
FNF{i}      = sum(confmatrixF{i},2);
FNF{i}      = FNF{i}(:)-TPF{i}(:);
TNF{i}      = sum(sum(confmatrixF{i}))-(TPF{i}+FPF{i}+FNF{i});
accF(i) = sum(TPF{i}./(TPF{i}+FPF{i}+TNF{i}+FNF{i}));
RIF(i) = sum((TPF{i}+TNF{i})./(TPF{i}+FPF{i}+TNF{i}+FNF{i}))/c;
PrecisionF{i} = TPF{i}./(TPF{i}+FPF{i});
RecallF{i} = TPF{i}./(TPF{i}+FNF{i});
F1Fm     = (2*PrecisionF{i}.*RecallF{i})./(PrecisionF{i}+RecallF{i});
F1Fm(isnan(F1Fm)) = 0;
F1F(i)     = sum(F1Fm)/c;
MIhatF(i) = MutualInfo(gnd,resF{i});
TPFs      = sum(TPF{i})/c;
FPFs      = sum(FPF{i})/c;
TNFs      = sum(TNF{i})/c;
FNFs      = sum(FNF{i})/c;
PrecisionFs(i)= sum(PrecisionF{i})/c;
RecallFs(i)   = sum(RecallF{i})/c;
end
 

% fip=  fopen(['results/LapL21selfrepresentation-ORL-',date,'.txt'],  'a+');
%     
%        
%     fprintf( fip, ['-------------------------------------LapL21selfrepresentation-ORL- dataset ----------------------------- \n']);
% fprintf( fip,'%s,lam1= %d,lam2= %d,L21accM=%.2f,L21accstd=%.2f,L21Mi=%.2f,L21Mistd=%.2f\n  ',datestr(now,31),lam1,lam2,mean(acc2)*100,std(acc2)*100,mean(MIhat2)*100,std(MIhat2)*100); 
% fprintf(  fip,'%s,lam1= %d,lam2= %d,F_L21accM=%.2f,F_L21accstd=%.2f,F_L21Mi=%.2f,F_L21Mistd=%.2f\n  ',datestr(now,31),lam1,lam2,mean(accF)*100,std(accF)*100,mean(MIhatF)*100,std(MIhatF)*100); 
% 
% fclose(fip);
cellnames = ['A',num2str(k),':V',num2str(k)];
result = [lam1,lam2,mean(TPs),std(TPs),mean(TNs),std(TNs),mean(FPs),std(FPs),mean(FNs),std(FNs),mean(acc)*100,std(acc)*100,...
                             mean(MIhat)*100,std(MIhat)*100,mean(Precisions),std(Precisions),mean(Recalls),std(Recalls),...
                             mean(F1),std(F1),mean(RI),std(RI)];
resultF= [lam1,lam2,mean(TPFs),std(TPFs),mean(TNFs),std(TNFs),mean(FPFs),std(FPFs),mean(FNFs),std(FNFs),mean(accF)*100,std(accF)*100,...
                             mean(MIhatF)*100,std(MIhatF)*100,mean(PrecisionFs),std(PrecisionFs),mean(RecallFs),std(RecallFs),...
                             mean(F1F),std(F1F),mean(RIF),std(RIF)];
xlswrite('AR_SGSSC', result,1,cellnames)
%k=k+1;
xlswrite('AR_SGSSC', resultF,2,cellnames)
k=k+1;
    end
end
% B = (abs(B)+ abs(B'))./2; 
% [acc MIhat] = computeACCNMI(B,gnd',rep);