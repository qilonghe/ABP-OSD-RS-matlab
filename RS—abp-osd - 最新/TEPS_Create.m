function T = TEPS_Create(k,Osd_i)
n = k;
TEPS = zeros(1,1); %初始化结果矩阵
idx = 1;
for Osd_i = 1
    iCns = nchoosek(1:n, Osd_i); %nchoosek 从n个元素中无顺序选出k个，穷尽所有可能情况
    for j = 1:size(iCns, 1)
        TEPS(idx, iCns(j, :)) = 1;
          idx = idx + 1;
    end
end
[m1,n1] = size(TEPS);
TEPS(m1+1,:) = zeros(1,k);
T =TEPS;
end
