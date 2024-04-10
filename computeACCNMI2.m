function [AC NMI L] = computeACCNMI2(B,gnd,rep)

c = length( unique( gnd ) ) ;
%Z = ( abs(T) + abs(T') ) / 2 ;
accAvg = [];
NMIAvg = [];
AC = [];
NMI =[];

for j = 1 : rep
% idx = clu_ncut(Z,nCluster) ;
% acc = compacc(idx,gnd); 
% accAvg = accAvg+acc;

L = BuildAdjacency(thrC(B),0); % 0.95
% --- perform Normalized Symmetric spectral clustering
rand('twister',5489);
P_label = SC(L,c);
label = reshape(P_label,1,[]);

res = bestMap(gnd,label);
acc = length(find(gnd == res))/length(gnd);
% accAvg = accAvg+acc;
accAvg = [accAvg acc];

MIhat = MutualInfo(gnd,res);
% NMIAvg = NMIAvg+ MIhat;
NMIAvg = [NMIAvg  MIhat];
end

AC.mean  = sum(accAvg)/rep;
AC.std = std(accAvg) ;

NMI.mean = sum(NMIAvg) /rep;
NMI.std = std(NMIAvg);

end