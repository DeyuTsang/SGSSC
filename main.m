close all; clc;clear
addpath('funs/')
addpath('results/')
addpath('data/')
load('ORL_32x32.mat')
% X   = X';
% gnd = X_label';
X   = double(X);
fea = NormalizeFea(X',0);
% options = [];
% options.ReducedDim = 800;
% [eigvector,eigvalue,meanData,new_data] = PCA(fea',options);
for lam1 = 1%[1e-3,2*1e-3,5*1e-3,8*1e-3,1e-2,2*1e-2,5*1e-2,8*1e-2,1e-1,0.2,0.5,0.8,1,2,5,8,10,50,100]
    for lam2 = 1
r=1;
c=max(gnd);
[Z2,F2] = LapdoubleL21SR(fea,lam1,lam2,c);
% [Z] = L21SR(fea,lambda,r);
[Z1] = doubleL21SR(fea,lam1);
% [Z2] = LRR(new_data',lambda);
% [Z,E] = lrra(X,X,lambda,1);
% figure;imshow(Z,[])
% figure;show((fea*Z)')
% figure;plot(Z(:,11),'b')


L1 = BuildAdjacency(thrC(Z1),1);
L2 = BuildAdjacency(thrC(Z2),1);
for i = 1:10
P_label1{i} = SC(L1,c);
label1{i} = reshape(P_label1{i},1,[]);
res1{i} = bestMap(gnd,label1{i});
acc1(i) = length(find(gnd == res1{i}))/length(gnd);
MIhat1(i) = MutualInfo(gnd,res1{i});


P_label2{i} = SC(L2,c);
label2{i} = reshape(P_label2{i},1,[]);
res2{i} = bestMap(gnd,label2{i});
acc2(i) = length(find(gnd == res2{i}))/length(gnd);
MIhat2(i) = MutualInfo(gnd,res2{i});

end
 

fip=  fopen(['results/L21selfrepresentation-ORL-PCA800-',date,'.txt'],  'a+');
    
       
    fprintf( fip, ['-------------------------------------L21selfrepresentation-ORL-PCA800- dataset ----------------------------- \n']);
fprintf( fip,'%s,lambda= %d,L21accM=%.2f,L21accstd=%.2f,L21Mi=%.2f,L21Mistd=%.2f\n  ',datestr(now,31),lam1,mean(acc1)*100,std(acc1)*100,mean(MIhat1)*100,std(MIhat1)*100); 
fprintf( fip,'%s,lambda= %d,Lap_accM=%.2f,Lap_accstd=%.2f,LarMi=%.2f,Lap_Mistd=%.2f\n  ',datestr(now,31),lam1,mean(acc2)*100,std(acc2)*100,mean(MIhat2)*100,std(MIhat2)*100); 

fclose(fip);
    end
end
% B = (abs(B)+ abs(B'))./2; 
% [acc MIhat] = computeACCNMI(B,gnd',rep);