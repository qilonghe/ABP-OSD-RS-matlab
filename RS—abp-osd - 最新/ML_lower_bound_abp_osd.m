function fer = ML_lower_bound_abp_osd(u_out,msg,receive)
fer = 0;
dist_dec2rec = sum(abs(u_out-receive));%译出码字和接收码字之间的距离
dist_dec2msg = sum(abs(u_out-msg));%译出码字和发送码字之间的距离
if  dist_dec2rec < dist_dec2msg %如果译出码字和接收码字之间的距离比译出码字和发送码字之间的距离小，那么即是是最大似然译码也会出错
   fer = 1;
end
end