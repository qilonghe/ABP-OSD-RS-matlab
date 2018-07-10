function llr = Create_Extllr(H3,llr,m,n,a_abp)
A = cell(1,n);         
B = cell(1,m);        

for vn=1:n
    A{vn} = find(H3(:,vn))'; % �����ڵ���У��ڵ������
end
for cn=1:m
    B{cn} = find(H3(cn,:));% У��ڵ��б����ڵ������
end
L = zeros(1,n);           % ���ÿ�������ڵ��llrֵ
M = ones(m,1) * llr;   % ��ʼ��llr��m*n
for cn=1:m
    for vn=B{cn}                   
        inM = M(cn,B{cn}(B{cn}~=vn)); % ���˱��������������ͬһ��У��ڵ������ӵı����ڵ�                  
        E(cn,vn) = 2*atanh_proceed(prod(tanh(inM/2))); %����������Ϣ               
    end            
end 
 % ����ÿһ�������ڵ㣬����Ϣ���ϱ�����Ϣ�õ��µ�llr
for vn=1:n            
    L(vn) = a_abp*sum(E(A{vn},vn))+llr(vn);             
end 
   llr = L;
end