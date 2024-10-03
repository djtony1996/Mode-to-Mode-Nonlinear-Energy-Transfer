clear,clc

a1 = load('M_590_52566252.mat');
a2 = load('M_590_52566254.mat');
a3 = load('M_590_52566298.mat');
a4 = load('M_590_52566299.mat');
a5 = load('M_590_52566300.mat');
a6 = load('M_590_52566302.mat');
a7 = load('M_590_52592832.mat');
a8 = load('M_590_52592960.mat');
a9 = load('M_590_52592962.mat');
a10= load('M_590_52592977.mat');
a11= load('M_590_52592991.mat');
a12= load('M_590_52593003.mat');
a13= load('M_590_52723928.mat');
%%
M_4d_avg = (a1.M_4d_avg + a2.M_4d_avg + a3.M_4d_avg + a4.M_4d_avg + a5.M_4d_avg + a6.M_4d_avg + a7.M_4d_avg + a8.M_4d_avg + a9.M_4d_avg + a10.M_4d_avg + a11.M_4d_avg + a12.M_4d_avg + a13.M_4d_avg) ./ 13;


%%

save('M_590_avg_final.mat','M_4d_avg','-v7.3')