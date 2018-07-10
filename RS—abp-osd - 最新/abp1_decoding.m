%% abp�㷨�����ڲ�ͬ����������ӵ�Ӱ��ܴ�Ҫ�ʵ�ѡȡ��
%% Ŀǰ���Թ��������rs���ǿ���Ӧ��
% hql
% 2017/12
function u_output = abp1_decoding(H,llr)
[m,n] = size(H);%m��У�鷽����Ŀn-k��nΪ�ܱ�����Ŀ
J = 1;
a_abp = 0.1;%abp������������
MaxIter_abp = 20;%abp����������
while(1)
[x,pos] = sort(abs(llr));%��llrֵ��С��������
H1 = H(:,pos);%��У�����Ҳһ������
%% �õ���H1���и�˹��Ԫ
H2 = bin_gauss_solve(H1);%�������ζȴ�С����У������˹��Ԫ
%% �ҳ�H2��ǰk���м��������޹�������λ
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
%% ��ǰk������������޹���ͺ�n-k��������ϳ����ɾ��������ֶ������ζȴ�С����
for i = 1:m
     H2_first(:,i) = H2(:,H2_pos_rk(i));
end
H2_pos1 = setdiff(pos,H2_pos_rk);
for i = 1:n-m
    H2_second(:,i) = H2(:,H2_pos1(i));
end
H2_total_pos = [H2_pos_rk,H2_pos1];
H2_u1 = [H2_first,H2_second];%�õ�ǰn-k���������޹����������Ұ������ζȴ�С�������򣬺�k��Ҳ�����ζȴ�С����--У�����ڶ�������
H3 = bin_gauss_solve(H2_u1);%�ٽ���һ�θ�˹��Ԫ
%��H3������˳��ԭ��ȥ
for i = 1:n
    H3_a(:,H2_total_pos(i)) = H3(:,i);
end
for i = 1:n
    H3_b(:,pos(i)) = H3_a(:,i);
end
%% �õ��µ�����Ϣ����llr
new_llr = Create_Extllr(H3_b,llr,m,n,a_abp);
llr = new_llr;

%% ��llr����Ӳ�о�
new_llr(new_llr<0)= 0;
new_llr(new_llr>0)= 1;
if(~any(mod(H*new_llr',2)))%�������ɹ�
    u_output = new_llr;
    return;
end
if (J == MaxIter_abp)
    u_output = new_llr;
    return;
end
J = J+1;%����ʧ�����µ�llr�ٴ��ݸ�osd������һ������
end
end