function out = expand_element(inH_element,B,m,order,Flagg)
out = zeros(m,m);
flag = 0;
for i = 1:m
    if inH_element == 0
        out = zeros(m,m);
        return
    else if mod(order(inH_element)+1+flag,2^m-1) ==0
        j = 2^m-1;
    else
        j = mod(order(inH_element)+1+flag,2^m-1);
        end
    end
    out(:,i)= B(:,j);
    flag = flag+1;
end
    if Flagg == 1 %%�������������G����Ҫ��չ��Ҫ��ʸ����չת�ã���
        out = out';
    end
end