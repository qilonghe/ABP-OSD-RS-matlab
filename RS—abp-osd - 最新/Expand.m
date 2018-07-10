function outH = Expand(inH,alpha,order,m,n,k)
B = zeros(m,n);
Flagg = 0;
[row col] = size(inH);
if row == n-k
    ExpandHorG = zeros((n-k)*m,n*m);%如果传进来的是H矩阵
else
    ExpandHorG = zeros(k*m,n*m);%如果是G矩阵
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