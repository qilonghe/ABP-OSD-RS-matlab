function [mult_ans]=multiply(a,b,m,alpha,order);
mult=0;
if  a*b==0
    mult_ans=0;
else
    a_order=order(a);%a的值对应的是order的列数，而该列中的元素即a对应的阶数
    b_order=order(b);
    mult_ans=mod(a_order+b_order,2^m-1);
    mult_ans=alpha(mult_ans+1);%mult的值对应的是alpha的列数，元素的阶数，而该列中的元素即为该元素的十进制值
end
%乘法结束时值表示 各阶数对应的十进制数值