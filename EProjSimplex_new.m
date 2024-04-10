function [x ft] = EProjSimplex_new(v, k)

%
%% Problem
%
%  min  1/2 || x - v||^2
%  s.t. x>=0, 1'x=1
%                   ������Ż�ʽ���൱��PCAN��p980 (27) or p979 (15)

if nargin < 2
    k = 1;
end;

ft=1;   % ��¼ѭ���˶��ٴ�
n = length(v);

%��ʵ�������v0��������s_ij�Ľ���ˣ������Ͽ���ֹ�������ⲽ��ֻ�����v0�������ŵģ�������v0�е�Զ������ĵ��������ٶ�������һЩ������Ӧ��������v0��������һ�����ŵ�ֵ������ʲô�ض����Ż��㷨
v0 = v-mean(v) + k/n;  %PCAN:p980 (30);����-mean(v)+k/n ��ʾ(30)ʽ���еĦǣ�����(30)ʽ�еĵ�һ����ʽ��v0Ӧ�õ���1�Ŷ�,v0��Ϊ����Ԫ�ؾ���s_ij.  1'x=1
%vmax = max(v0);
vmin = min(v0);
if vmin < 0
    f = 1;
    lambda_m = 0;
    while abs(f) > 10^-10
        v1 = v0 - lambda_m;
        posidx = v1 > 0;
        npos = sum(posidx);     %��ǰv0���ж�����Ч��s_ij
        g = -npos;
        f = sum(v1(posidx)) - k;        %��ǰ��1�Ĳ�ֵ
        lambda_m = lambda_m - f/g;      %���߲�ֵf ������lambda_m ���ٸ���lambda_m������v1,����v1�����ͱƽ�1����Խ��lambda_mԽ��,v0�ͼ�ȥԽ�������ͬʱs_ij�ĸ���g����̫��
        ft=ft+1;
        if ft > 100
            x = max(v1,0);
            break;
        end;
    end;
    x = max(v1,0);      %ֻȡ������Ĳ���Ϊs_ij

else
    x = v0;
end;