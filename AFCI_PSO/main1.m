clc; clear all;
func_num=3;
D=18;
Xmin=-4;
Xmax=4;
pop_size=10000;
iter_max=1500;
% group_num = 2;
group_num = 500;
runs=1;
fhd=str2func('cec19_func');

for i = 1:runs
    fileID=fopen('result.txt','a');
    [gbest_b,gbestval_b,g_index,Iter,FES]= SC_Func(fhd,D,pop_size,iter_max,group_num,Xmin,Xmax,func_num)
    % gbestval_b = gbestval_b*10^9;
    fprintf(fileID,'%d th epochs:%d iters, %d fes\n',i,Iter,FES);
    fprintf(fileID,'%f\n',gbestval_b);
    gbesteval(:,:,i) = gbestval_b;
    ite(:,:,i) = Iter;
    FE(:,:,i) = FES;
    fclose(fileID);
    % gbest(:,:,i) = gbest_b;
end

    s_g = squeeze(gbesteval);
    % tongji(1,1) = min(s_g,[],[1 2]);
    [tongji(1,1),index] = min(s_g);
    tongji(2,1) = max(s_g,[],[1 2]);
    tongji(3,1) = mean(s_g,[1 2]);
    tongji(4,1) = median(s_g,[1 2]);
    tongji(5,1) = var(s_g,0,[1 2]);

save res.txt -ascii tongji;
% for i=1:29
%     func_num=(i<2) * i + (i>=2) * (i+1);
%     for j=1:runs
%         i,j,
%         [gbest_b,gbestval_b,FES]= SC_Func(fhd,D,pop_size,iter_max,group_num,Xmin,Xmax,func_num);
% %         xbest(j,:)=gbest_b;
% %         fbest(i,j)=gbestval_b;
% %         fbest(i,j)
%     end
% %     f_mean(i)=mean(fbest(i,:));
% end
