%% osd������abp�����㷨(���rs����Ƽ���)
% hql
% 2017/12
function u_output = abp_osd_decoding(llr,H,G,r_abp_osd,u_output_test)
[k,n] = size(G);
MaxIter_abp = 20;%abp����������
a_abp = 0.1;%abp������������
J = 1; %abp��ʼ����������
aaa=0;
corre_num = zeros(1,MaxIter_abp);
u_output_osd = zeros(MaxIter_abp,n);
%% %%%%%%%%%%%%%%%%%%%  OSD����׶�  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% osd����
while(1)
    [u_output_u2,fin_d,correlate_val] =osd2_decoding(llr,G,r_abp_osd);
    u_output_osd(J,:) = u_output_u2;%osd�������
    %% �ж�osd�����Ƿ��ҵ����������֣�����ҵ���ֱ��������������ж��Ƿ񵽴�abp������������������������abp����
    corre_num(J) = correlate_val;
    if  fin_d == 1 % �ҵ����������֡�
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
    %% %%%%%%%%%%%%%%%%%%% abp����׶� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%������ֱ�ӵ���abp�㷨����Ϊ������abpֻ����һ�ξͽ�llr����osd���������
    [mu,nu] = size(H);
    [x,pos] = sort(abs(llr));%��llrֵ��С��������
    H1 = H(:,pos);%��У�����Ҳһ������
    %% �õ���H1���и�˹��Ԫ
    H2 = bin_gauss_solve(H1);%�������ζȴ�С����У������˹��Ԫ
    %% �ҳ�H2��ǰk���м��������޹�������λ
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
    %% ��ǰk������������޹���ͺ�n-k��������ϳ����ɾ��������ֶ������ζȴ�С����
    for i = 1:mu
         H2_first(:,i) = H2(:,H2_pos_rk(i));
    end
    H2_pos1 = setdiff(pos,H2_pos_rk);
    for i = 1:nu-mu
        H2_second(:,i) = H2(:,H2_pos1(i));
    end
    H2_total_pos = [H2_pos_rk,H2_pos1];
    H2_u1 = [H2_first,H2_second];%�õ�ǰn-k���������޹����������Ұ������ζȴ�С�������򣬺�k��Ҳ�����ζȴ�С����--У�����ڶ�������
    H3 = bin_gauss_solve(H2_u1);%�ٽ���һ�θ�˹��Ԫ
    %��H3������˳��ԭ��ȥ
    for i = 1:nu
        H3_a(:,H2_total_pos(i)) = H3(:,i);
    end
    for i = 1:nu
        H3_b(:,pos(i)) = H3_a(:,i);
    end
    %% �õ��µ�����Ϣ����llr
    new_llr = Create_Extllr(H3_b,llr,mu,nu,a_abp);
    llr = new_llr;
    %% ���µ�llr����Ӳ�о�
%    new_llr(new_llr<0)= 0;
%    new_llr(new_llr>0)= 1;
%    if(~any(mod(H*new_llr',2)))%�������ɹ�
%        u_output = new_llr;
%        if sum(u_output~=u_output_test)
%            aaa=1;
%            u_output = u_output_test;%���abp����Ľ����osd��ͬ����ôΪ��ֹabp�����˺Ϸ�����ȴ������ȷ���֣��ͻ������osd�������
%            fprintf('aaa=',aaa);
%        end
%        return;
%    end
    J = J+1;%����ʧ�����µ�llr�ٴ��ݸ�osd������һ������
end


end