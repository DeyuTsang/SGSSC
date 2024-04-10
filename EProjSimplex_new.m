function [x ft] = EProjSimplex_new(v, k)

%
%% Problem
%
%  min  1/2 || x - v||^2
%  s.t. x>=0, 1'x=1
%                   上面的优化式子相当于PCAN：p980 (27) or p979 (15)

if nargin < 2
    k = 1;
end;

ft=1;   % 记录循环了多少次
n = length(v);

%其实下面这个v0就是求解的s_ij的结果了，理论上可以止于下面这步。只是这个v0不是最优的，甚至离v0有点远，后面的迭代就是再对它进行一些调整，应该是想在v0附近搜索一个更优的值。不是什么特定的优化算法
v0 = v-mean(v) + k/n;  %PCAN:p980 (30);后面-mean(v)+k/n 表示(30)式的中的η，根据(30)式中的第一个等式，v0应该等于1才对,v0中为正的元素就是s_ij.  1'x=1
%vmax = max(v0);
vmin = min(v0);
if vmin < 0
    f = 1;
    lambda_m = 0;
    while abs(f) > 10^-10
        v1 = v0 - lambda_m;
        posidx = v1 > 0;
        npos = sum(posidx);     %当前v0中有多少有效的s_ij
        g = -npos;
        f = sum(v1(posidx)) - k;        %当前与1的差值
        lambda_m = lambda_m - f/g;      %更具差值f 来调整lambda_m ，再更具lambda_m来调整v1,想让v1正数和逼近1。差越大，lambda_m越大,v0就减去越大的数。同时s_ij的个数g不能太少
        ft=ft+1;
        if ft > 100
            x = max(v1,0);
            break;
        end;
    end;
    x = max(v1,0);      %只取大于零的部分为s_ij

else
    x = v0;
end;