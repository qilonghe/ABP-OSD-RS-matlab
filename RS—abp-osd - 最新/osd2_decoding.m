%%osd算法
%%hql-2017/12
function [u_output,fin_d,correlate_val] = osd2_decoding(llr,G,r_osd)%fin_d==0即没找到最大概似码字，输出的是相关性最大码字，为1为找到了最大概似码字
[k,n] = size(G);
posSot = 1:n;
%%%%%%%%%%%%%%%%%%%  OSD译码  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% osd信任度处理
[xx,posSot] = sort(abs(llr),'descend');%将llr值从大到小排序
llrSort = llr(posSot);%得到排序后的llr值，信任度从大到小排序

r1 = r_osd(:,posSot);%将接收值r按照信任度排序--r第一次排序。
G1 = G(:,posSot);%将生成矩阵也一样排序--生成矩阵第一次排序


%% 得到的G2进行高斯消元（只进行行变换）

G2= bin_gauss_solve(G1);%根据信任度从大到小将生成矩阵高斯消元
%% 找出G2中前k个列极大线性无关向量组位
p = 0;
G2_pos = [];
for i = 1:k
    j = i;
      while( 1 )
        if(G2(i,j)==1)
           G2_pos(i) = j;
           break;
        end 
        j = j+1;
        while(find(G2_pos == j) ~= 0)
            j = j+1;
        end
      end     
end
G2_pos_rk = sort(G2_pos);
%% 将前k个列最大线性无关组和后n-k个重新组合成生成矩阵（两部分都是信任度从大到小）
for i = 1:k
     G2_first(:,i) = G2(:,G2_pos_rk(i));
end
G2_pos1 = setdiff(posSot,G2_pos_rk);
for i = 1:n-k
    G2_second(:,i) = G2(:,G2_pos1(i));
end
G2_total_pos = [G2_pos_rk,G2_pos1];
G2_u1 = [G2_first,G2_second];%得到前k个是线性无关列向量，且按照信任度从大到小排序，后n-k个也是信任度从大到小--生成矩阵第二次排序
G3 = bin_gauss_solve(G2_u1);%再进行一次高斯消元
r2 = r1(:,G2_total_pos);%接收向量r同样根据生成矩阵排序--r第二次排序



%% 得到前k个信任度最大的比特，并且进行硬判决
MRIPPos = r2(1:k);%得到MRIP
MRIPPos(MRIPPos<0) = 0; %将MRIP硬判决
MRIPPos(MRIPPos>0) = 1;

%% 选择osd译码阶数并得到所有错误图案
Osd_i = 1; %选择osd阶数,目前1是秒出的，如果2阶，当k较大时，复杂度较高，差别不大，一般用1阶或者0阶。
TEPS = TEPS_Create(k,Osd_i);
%得到修正后所有码字组合
for i = 1:k+1
    MRIP(i,:) = xor(TEPS(i,:),MRIPPos);
end
%将接收值r进行硬判决
r3 = zeros(1,n);
r3(r2<0)= 0;
r3(r2>0)= 1;

%% 根据前k个信任度最大的比特重新生成所有码字
c_possible_total = mod(MRIP*G3, 2);%得到所有可能的码字集合

%% 在候选码字中进行寻找到正确码字
Dmin = 7;%针对于不同码，有不同的最小码距，rs码的（7,3,5），（31,25,7）,(15,7,9)
fin_d = 0;
for i = 1:k+1
    D1 = find(c_possible_total(i,:) ~= r3);%码字中，与接收信号不同的位置集合；
    D0 = setdiff([1:n],D1);%码字中，与接收信号相同的位置结合
    D0_size= size(D0,2);
    D1_size = size(D1,2);%码字中不同于接收信号的个数
   % if  Dmin > D1_size %只有小于dmin的不同码字个数的时候才有必要来判断是否为最大概似码字
        Dist_Code(i) = sum(abs(r2(find(c_possible_total(i,:) ~= r3)))); %得到接收信号和码字间的相关性差异λ
        D0_A = D0(D0_size-(Dmin-D1_size)+1:D0_size);%得到D0中信任度最小的dmin-D0size的位置集合
        G = sum(abs(r2(D0_A)));%定义的比较系数
        if Dist_Code(i) < G %如果成立，就为最大概似码字，直接输出码字
        c_possible = c_possible_total(i,:);
        fin_d = 1;
        correlate_val = 1000;%如果找到最大概似码字了，把相关性置为最大
        break;
       end
 %   end   
    for ii = 1:n
       theta(ii) = r2(ii)*(2*c_possible_total(i,ii)-1);
    end
    theta_A(i) = sum(theta);%得到当前码字的相关性theta

end
if fin_d == 0                               
   [val,pos] = sort(theta_A,'descend');%如果找不到最大概似码字，就输出相关性最大的
   c_possible = c_possible_total(pos(1),:);
   correlate_val = val(1);%如果没有找到最大概似码字就保留相关性最大的值以及对应的码字，
end
%% 将得到的最小距离的码字顺序还原回去
for i = 1 : n
   u_output_u1(:,G2_total_pos(i)) = c_possible(i);%两个for不能进行合并
end
for i = 1 : n
    u_output_u2(:,posSot(i)) = u_output_u1(i);
end
u_output = u_output_u2;%

end