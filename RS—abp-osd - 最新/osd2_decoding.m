%%osd�㷨
%%hql-2017/12
function [u_output,fin_d,correlate_val] = osd2_decoding(llr,G,r_osd)%fin_d==0��û�ҵ����������֣�������������������֣�Ϊ1Ϊ�ҵ�������������
[k,n] = size(G);
posSot = 1:n;
%%%%%%%%%%%%%%%%%%%  OSD����  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% osd���ζȴ���
[xx,posSot] = sort(abs(llr),'descend');%��llrֵ�Ӵ�С����
llrSort = llr(posSot);%�õ�������llrֵ�����ζȴӴ�С����

r1 = r_osd(:,posSot);%������ֵr�������ζ�����--r��һ������
G1 = G(:,posSot);%�����ɾ���Ҳһ������--���ɾ����һ������


%% �õ���G2���и�˹��Ԫ��ֻ�����б任��

G2= bin_gauss_solve(G1);%�������ζȴӴ�С�����ɾ����˹��Ԫ
%% �ҳ�G2��ǰk���м��������޹�������λ
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
%% ��ǰk������������޹���ͺ�n-k��������ϳ����ɾ��������ֶ������ζȴӴ�С��
for i = 1:k
     G2_first(:,i) = G2(:,G2_pos_rk(i));
end
G2_pos1 = setdiff(posSot,G2_pos_rk);
for i = 1:n-k
    G2_second(:,i) = G2(:,G2_pos1(i));
end
G2_total_pos = [G2_pos_rk,G2_pos1];
G2_u1 = [G2_first,G2_second];%�õ�ǰk���������޹����������Ұ������ζȴӴ�С���򣬺�n-k��Ҳ�����ζȴӴ�С--���ɾ���ڶ�������
G3 = bin_gauss_solve(G2_u1);%�ٽ���һ�θ�˹��Ԫ
r2 = r1(:,G2_total_pos);%��������rͬ���������ɾ�������--r�ڶ�������



%% �õ�ǰk�����ζ����ı��أ����ҽ���Ӳ�о�
MRIPPos = r2(1:k);%�õ�MRIP
MRIPPos(MRIPPos<0) = 0; %��MRIPӲ�о�
MRIPPos(MRIPPos>0) = 1;

%% ѡ��osd����������õ����д���ͼ��
Osd_i = 1; %ѡ��osd����,Ŀǰ1������ģ����2�ף���k�ϴ�ʱ�����ӶȽϸߣ���𲻴�һ����1�׻���0�ס�
TEPS = TEPS_Create(k,Osd_i);
%�õ������������������
for i = 1:k+1
    MRIP(i,:) = xor(TEPS(i,:),MRIPPos);
end
%������ֵr����Ӳ�о�
r3 = zeros(1,n);
r3(r2<0)= 0;
r3(r2>0)= 1;

%% ����ǰk�����ζ����ı�������������������
c_possible_total = mod(MRIP*G3, 2);%�õ����п��ܵ����ּ���

%% �ں�ѡ�����н���Ѱ�ҵ���ȷ����
Dmin = 7;%����ڲ�ͬ�룬�в�ͬ����С��࣬rs��ģ�7,3,5������31,25,7��,(15,7,9)
fin_d = 0;
for i = 1:k+1
    D1 = find(c_possible_total(i,:) ~= r3);%�����У�������źŲ�ͬ��λ�ü��ϣ�
    D0 = setdiff([1:n],D1);%�����У�������ź���ͬ��λ�ý��
    D0_size= size(D0,2);
    D1_size = size(D1,2);%�����в�ͬ�ڽ����źŵĸ���
   % if  Dmin > D1_size %ֻ��С��dmin�Ĳ�ͬ���ָ�����ʱ����б�Ҫ���ж��Ƿ�Ϊ����������
        Dist_Code(i) = sum(abs(r2(find(c_possible_total(i,:) ~= r3)))); %�õ������źź����ּ������Բ����
        D0_A = D0(D0_size-(Dmin-D1_size)+1:D0_size);%�õ�D0�����ζ���С��dmin-D0size��λ�ü���
        G = sum(abs(r2(D0_A)));%����ıȽ�ϵ��
        if Dist_Code(i) < G %�����������Ϊ���������֣�ֱ���������
        c_possible = c_possible_total(i,:);
        fin_d = 1;
        correlate_val = 1000;%����ҵ������������ˣ����������Ϊ���
        break;
       end
 %   end   
    for ii = 1:n
       theta(ii) = r2(ii)*(2*c_possible_total(i,ii)-1);
    end
    theta_A(i) = sum(theta);%�õ���ǰ���ֵ������theta

end
if fin_d == 0                               
   [val,pos] = sort(theta_A,'descend');%����Ҳ������������֣���������������
   c_possible = c_possible_total(pos(1),:);
   correlate_val = val(1);%���û���ҵ����������־ͱ������������ֵ�Լ���Ӧ�����֣�
end
%% ���õ�����С���������˳��ԭ��ȥ
for i = 1 : n
   u_output_u1(:,G2_total_pos(i)) = c_possible(i);%����for���ܽ��кϲ�
end
for i = 1 : n
    u_output_u2(:,posSot(i)) = u_output_u1(i);
end
u_output = u_output_u2;%

end