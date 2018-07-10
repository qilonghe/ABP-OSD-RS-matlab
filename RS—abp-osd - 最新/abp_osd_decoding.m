%% osd译码与abp级联算法(针对rs码设计级联)
% hql
% 2017/12
function u_output = abp_osd_decoding(llr,H,G,r_abp_osd,u_output_test)
[k,n] = size(G);
MaxIter_abp = 20;%abp迭代最大次数
a_abp = 0.1;%abp迭代抑制因子
J = 1; %abp初始化迭代次数
aaa=0;
corre_num = zeros(1,MaxIter_abp);
u_output_osd = zeros(MaxIter_abp,n);
%% %%%%%%%%%%%%%%%%%%%  OSD译码阶段  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% osd译码
while(1)
    [u_output_u2,fin_d,correlate_val] =osd2_decoding(llr,G,r_abp_osd);
    u_output_osd(J,:) = u_output_u2;%osd译码输出
    %% 判断osd译码是否找到最大概似码字，如果找到则直接译码结束，再判断是否到达abp译码的最大迭代次数。否则进行abp译码
    corre_num(J) = correlate_val;
    if  fin_d == 1 % 找到最大概似码字。
        u_output = u_output_osd(J,:);
        if sum(u_output~=u_output_test)
            aaa=1;
            fprintf('aaa=',aaa);
        end
        return;
    end
    if  J == MaxIter_abp
        [val_l pos_s] = sort(corre_num,'descend');

        Max_Correlate_num = pos_s(1);
        u_output = u_output_osd(Max_Correlate_num,:);
        if sum(u_output~=u_output_test)
            aaa=1;
            fprintf('aaa=',aaa);
        end
        return;
    end
    %% %%%%%%%%%%%%%%%%%%% abp译码阶段 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%（不能直接调用abp算法，因为级联中abp只迭代一次就将llr给了osd或者输出）
    [mu,nu] = size(H);
    [x,pos] = sort(abs(llr));%将llr值从小到大排序
    H1 = H(:,pos);%将校验矩阵也一样排序
    %% 得到的H1进行高斯消元
    H2 = bin_gauss_solve(H1);%根据信任度从小到大将校验矩阵高斯消元
    %% 找出H2中前k个列极大线性无关向量组位
    H2_pos = [];
    for i = 1:mu
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
    for i = 1:mu
         H2_first(:,i) = H2(:,H2_pos_rk(i));
    end
    H2_pos1 = setdiff(pos,H2_pos_rk);
    for i = 1:nu-mu
        H2_second(:,i) = H2(:,H2_pos1(i));
    end
    H2_total_pos = [H2_pos_rk,H2_pos1];
    H2_u1 = [H2_first,H2_second];%得到前n-k个是线性无关列向量，且按照信任度从小到大排序，后k个也是信任度从小到大--校验矩阵第二次排序
    H3 = bin_gauss_solve(H2_u1);%再进行一次高斯消元
    %将H3矩阵列顺序还原回去
    for i = 1:nu
        H3_a(:,H2_total_pos(i)) = H3(:,i);
    end
    for i = 1:nu
        H3_b(:,pos(i)) = H3_a(:,i);
    end
    %% 得到新的外信息更新llr
    new_llr = Create_Extllr(H3_b,llr,mu,nu,a_abp);
    llr = new_llr;
    %% 将新的llr进行硬判决
%    new_llr(new_llr<0)= 0;
%    new_llr(new_llr>0)= 1;
%    if(~any(mod(H*new_llr',2)))%如果译码成功
%        u_output = new_llr;
%        if sum(u_output~=u_output_test)
%            aaa=1;
%            u_output = u_output_test;%如果abp译码的结果和osd不同，那么为防止abp出现了合法码字却不是正确码字，就还是输出osd译的码字
%            fprintf('aaa=',aaa);
%        end
%        return;
%    end
    J = J+1;%译码失败则将新的llr再传递给osd继续下一次译码
end


end