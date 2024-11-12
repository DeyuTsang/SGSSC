function [TPs,FPs,TNs,FNs,RI,acc,Precisions,Recalls,F1,MIhat] = computeRes(gnd,res)
%UNTITLED2 此处提供此函数的摘要
%   此处提供详细说明
c = max(res);
confmatrix = cfmatrix(gnd, res, 1:c, 0);
TP      = diag(confmatrix);
FP      = sum(confmatrix);
FP      = FP(:)-TP(:);
FN      = sum(confmatrix,2);
FN      = FN(:)-TP(:);
TN      = sum(sum(confmatrix))-(TP+FP+FN);
RI = sum((TP+TN)./(TP+FP+TN+FN))/c;
acc = sum(TP./(TP+FP+TN+FN));
Precision = TP./((TP+FP)+eps);
Recall = TP./(TP+FN);
F1m     = (2*Precision.*Recall)./(Precision+Recall);
F1m(isnan(F1m)) = 0;
F1   = sum(F1m)/c;
MIhat = MutualInfo(gnd,res);
TPs      = sum(TP)/c;
FPs      = sum(FP)/c;
TNs      = sum(TN)/c;
FNs      = sum(FN)/c;
Precisions= sum(Precision)/c;
Recalls   = sum(Recall)/c;
end