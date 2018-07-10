%% abp算法，对于不同的码迭代因子的影响很大，要适当选取。
%% 目前尝试过汉明码和rs码是可以应用
% hql
% 2017/12
function u_output = abp1_decoding(H,llr)
[m,n] = size(H);%m是校验方程数目n-k。n为总比特数目
J = 1;
a_abp = 0.1;%abp迭代抑制因子
MaxIter_abp = 20;%abp迭代最大次数
while(1)
[x,pos] = sort(abs(llr));%将llr值从小到大排序
H1 = H(:,pos);%将校验矩阵也一样排序
%% 得到的H1进行高斯消元
H2 = bin_gauss_solve(H1);%根据信任度从小到大将校验矩阵高斯消元
%% 找出H2中前k个列极大线性无关向量组位
H2_pos = [];
for i = 1:m
    j = i;
      while( 1 )
        if(H2(i,j)==1)
           H2_pos(i) = j;
           break;
        end 
        j = j+1;
        while(find(H2_pos == j) ~= 0)
            j = j+1;
        end
      end     
end
H2_pos_rk = sort(H2_pos);
%% 将前k个列最大线性无关组和后n-k个重新组合成生成矩阵（两部分都是信任度从小到大）
for i = 1:m
     H2_first(:,i) = H2(:,H2_pos_rk(i));
end
H2_pos1 = setdiff(pos,H2_pos_rk);
for i = 1:n-m
    H2_second(:,i) = H2(:,H2_pos1(i));
end
H2_total_pos = [H2_pos_rk,H2_pos1];
H2_u1 = [H2_first,H2_second];%得到前n-k个是线性无关列向量，且按照信任度从小到大排序，后k个也是信任度从小到大--校验矩阵第二次排序
H3 = bin_gauss_solve(H2_u1);%再进行一次高斯消元
%将H3矩阵列顺序还原回去
for i = 1:n
    H3_a(:,H2_total_pos(i)) = H3(:,i);
end
for i = 1:n
    H3_b(:,pos(i)) = H3_a(:,i);
end
%% 得到新的外信息更新llr
new_llr = Create_Extllr(H3_b,llr,m,n,a_abp);
llr = new_llr;

%% 将llr进行硬判决
new_llr(new_llr<0)= 0;
new_llr(new_llr>0)= 1;
if(~any(mod(H*new_llr',2)))%如果译码成功
    u_output = new_llr;
    return;
end
if (J == MaxIter_abp)
    u_output = new_llr;
    return;
end
J = J+1;%译码失败则将新的llr再传递给osd继续下一次译码
end
end