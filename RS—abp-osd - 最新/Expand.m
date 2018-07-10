function outH = Expand(inH,alpha,order,m,n,k)
B = zeros(m,n);
Flagg = 0;
[row col] = size(inH);
if row == n-k
    ExpandHorG = zeros((n-k)*m,n*m);%�������������H����
else
    ExpandHorG = zeros(k*m,n*m);%�����G����
    Flagg = 1;
end
for i = 1:n
    B(:,i) = de2bi(alpha(i),m)';
end
if Flagg == 1
    x = k;
else
    x = n-k;
end
for i = 1:x
    for j = 1:n  
       ExpandHorG((i-1)*m+1:(i-1)*m+m,(j-1)*m+1:(j-1)*m+m)= expand_element(inH(i,j),B,m,order,Flagg);
    end
end
outH = ExpandHorG;
end