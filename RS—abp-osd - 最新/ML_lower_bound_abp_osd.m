function fer = ML_lower_bound_abp_osd(u_out,msg,receive)
fer = 0;
dist_dec2rec = sum(abs(u_out-receive));%������ֺͽ�������֮��ľ���
dist_dec2msg = sum(abs(u_out-msg));%������ֺͷ�������֮��ľ���
if  dist_dec2rec < dist_dec2msg %���������ֺͽ�������֮��ľ����������ֺͷ�������֮��ľ���С����ô�����������Ȼ����Ҳ�����
   fer = 1;
end
end