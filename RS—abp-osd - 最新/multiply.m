function [mult_ans]=multiply(a,b,m,alpha,order);
mult=0;
if  a*b==0
    mult_ans=0;
else
    a_order=order(a);%a��ֵ��Ӧ����order���������������е�Ԫ�ؼ�a��Ӧ�Ľ���
    b_order=order(b);
    mult_ans=mod(a_order+b_order,2^m-1);
    mult_ans=alpha(mult_ans+1);%mult��ֵ��Ӧ����alpha��������Ԫ�صĽ������������е�Ԫ�ؼ�Ϊ��Ԫ�ص�ʮ����ֵ
end
%�˷�����ʱֵ��ʾ ��������Ӧ��ʮ������ֵ