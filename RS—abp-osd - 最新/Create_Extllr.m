function llr = Create_Extllr(H3,llr,m,n,a_abp)
A = cell(1,n);         
B = cell(1,m);        

for vn=1:n
    A{vn} = find(H3(:,vn))'; % 变量节点中校验节点的坐标
end
for cn=1:m
    B{cn} = find(H3(cn,:));% 校验节点中变量节点的坐标
end
L = zeros(1,n);           % 存放每个变量节点的llr值
M = ones(m,1) * llr;   % 初始化llr，m*n
for cn=1:m
    for vn=B{cn}                   
        inM = M(cn,B{cn}(B{cn}~=vn)); % 除了本身以外的其他在同一个校验节点下连接的变量节点                  
        E(cn,vn) = 2*atanh_proceed(prod(tanh(inM/2))); %产生的外信息               
    end            
end 
 % 对于每一个变量节点，外信息加上本身信息得到新的llr
for vn=1:n            
    L(vn) = a_abp*sum(E(A{vn},vn))+llr(vn);             
end 
   llr = L;
end