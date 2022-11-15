function [gbest_b,gbesteval_b,g_index,Iter,fitcount] = SC_Func(fhd,Dimension,Particle_Number,Max_Gen,Group_num,VRmin,VRmax,varargin)

rand('state',sum(100*clock));
D = Dimension; Pop_size = Particle_Number; G_n = Group_num; Iter_times = Max_Gen;
inertia = 0.7;to_leader = 0.2;to_otherg = 0.1; threshold = 0.1;
C_hop = 0.1; 
C_Visit = 20; 
C_Recon = 40;

% C_hop = 0.0; 
% C_Visit= 10000000; 
% C_Recon = 10000000;
W1 = 0.8; W2 = 0.2; C_DE = 1; func_num = varargin{:};
Vmin = 0.5 * VRmin;
Vmax = 0.5 * VRmax;
F = 0.6;
CR = 0.8; stopeval = 1.000000000;
pos = bsxfun(@times, VRmin+(VRmax-VRmin), rand(Pop_size,D));
% pos=VRmin+(VRmax-VRmin).*rand(Pop_size,D);  % 初始的输入矩阵；生成的范围肯定是在[-100,100]之间，相当于原版算法中的X，Y定义
e = 5000 .* rand(1,Pop_size);
% e = feval(fhd,pos',varargin{:});          % 使用{:}补位 varargin，防止因为维度不匹配的原因，feval函数执行失败，且这里输出的是针对初始的输入矩阵的

fitcount = 0;
% fitcount = Pop_size;                        
vel = Vmin+2.*Vmax.*rand(Pop_size,D);       
%% 科学合作算法初始化
modd = mod(Pop_size, G_n);
if modd ~= 0 
    Pop_size = Pop_size + G_n - modd;
    pos = [pos;NaN((G_n-modd),D)];        
    e = [e,NaN(1,(G_n-modd))];            
    vel = [vel;NaN((G_n-modd),D)];
end
G_p_n = Pop_size/G_n;                       % 初始分堆每个种群中的个体的数量，考虑若不能整除的情况，判断有无余数；用全0补位
balance1 = ceil(G_p_n*2/3);                         % 小组平衡参数（应该是用于个体迁移所使用的）

% 分页：
g_pbest = permute(reshape(pos',D,G_p_n,G_n),[2 1 3]);                                   % 将pos矩阵分页
g_pvel = permute(reshape(vel',D,G_p_n,G_n),[2 1 3]);                                    % 将速度矩阵分页
g_pbesteval = reshape(e',G_p_n,1,G_n);                                                  % 这里就已经将适应度行向量转置成列向量
g_VRmin = repmat(VRmin, 1, 1, G_n);
g_VRmax = repmat(VRmax, 1, 1, G_n);
g_Vmin = repmat(Vmin, 1, 1, G_n);
g_Vmax = repmat(Vmax, 1, 1, G_n);

[gbesteval, gbestid] = min(g_pbesteval);                                                % 
for i = 1:G_n
    gbest(:,:,i) = g_pbest(gbestid(1,:,i), :,i);         % ok, 我找到了每页中的leader的行向量
    gbestrep(:,:,i) = repmat(gbest(:,:,i),G_p_n,1);      % 扩展leader行向量成为一个G_p_n行的矩阵
end

Vg = zeros(1,D,G_n);
Vgxy = zeros(1,D,G_n);                                                      
%% 迭代开始
   
% CMA_ES parameter setting
xstart = rand(D,1);
insigma = 0.3;
inopts.MaxIter = 1;
inopts.LBounds = VRmin;
inopts.UBounds = VRmax;
inopts.DispFinal = 0;
inopts.DispModulo = 0;
inopts.SaveVariables = 0;
inopts.SaveFilename = 0;
inopts.LogModulo = 0;
inopts.LogTime = 0;
inopts.LogFilenamePrefix = 0;
inopts.LogPlot = 0;

% arx = zeros(G_p_n,D,G_n);
% 种群融合加入进来
% 差分进化交叉操作
for j1 = 1:Iter_times
    ttmp = squeeze(g_pbesteval);
    M1 = min(ttmp(:));
    % M1 = min(g_pbesteval,[],[1 3]);
    
    % this is SRGO

    maxval = M1;                                                            
    j2 = 1:G_n;
    C_Group = to_otherg .* bsxfun(@rdivide, (maxval./gbesteval(:,:,j2)), (1+(maxval./gbesteval(:,:,j2))));
    %C_Group = to_otherg .* ((maxval./gbesteval(:,:,j2))./(1+(maxval./gbesteval(:,:,j2))));
    % 向下开始更新种群中个体的位置
    to_leader_t = 1 - to_leader - C_Group;
    g_pbest_tmp = g_pbest;  
    g_pvel_tmp = g_pvel;
    g_pvel = bsxfun(@times, inertia, g_pvel) + bsxfun(@times ,to_leader_t, bsxfun(@minus, gbest, g_pbest)) + bsxfun(@times, C_Group, Vg);
    % 更新速度和坐标
    g_pvel(bsxfun(@gt,g_pvel, g_Vmax)) = Vmax;
    g_pvel(bsxfun(@lt,g_pvel, g_Vmin)) = Vmin;
    
    % 更新坐标
    g_pbest = g_pbest + g_pvel;

    g_pbest = (g_pbest > g_VRmax) .* (g_VRmax - 0.25 .* (g_VRmax - g_VRmin) .* rand(size(g_pbest,1),D,G_n)) + (g_pbest < g_VRmin) .* (g_VRmin + 0.25 .* (g_VRmax - g_VRmin) .* rand(size(g_pbest,1),D,G_n))...
        + ((g_pbest >= g_VRmin) & (g_pbest <= g_VRmax)) .* g_pbest;
    % [g_pbest,g_pbesteval,g_pvel,fitcount] = select(g_pbest_tmp,g_pvel_tmp,g_pbest,g_pbesteval,g_pvel,G_n,func_num,fitcount);
    
    for j3 = 1:size(g_pbest,3)
        tmp = g_pbest(:,:,j3);
        tmp = rmmissing(tmp);
        tmp1 = feval(fhd, tmp', varargin{:});
        g_pbesteval(:,:,j3) = [tmp1'; NaN(size(g_pbesteval,1) - size(tmp1,2), 1)];
        fitcount = fitcount + size(tmp, 1);
    end
    
    if mod(j1, C_DE) == 0 && Iter_times - j1 > 10
        % g_dv = VRmin + (VRmax-VRmin) .* zeros(size(g_pvel));
        [g_pbest, g_pbesteval, g_pvel,fitcount] = DE_M_C_S(g_pbest,g_pbesteval,g_pvel,G_n,F,VRmax,VRmin,func_num,fitcount);
    end
    

      
    % Start hopping
    if G_n ~= 1
       
        [g_pbest, g_pbesteval, g_pvel] = hopping(D, g_pbest, g_pbesteval, g_pvel, gbestid, G_n, C_hop, balance1);
        
        % 开始更新新的leader
        [gbesteval, gbestid] = min(g_pbesteval);
        
        F_S = 0;
        F_V = zeros(1,G_n + 1);
        for i = 1:G_n
            gbest(:,:,i) = g_pbest(gbestid(1,:,i), :,i);         % ok, 我找到了每页中的leader的行向量
            gbestrep(:,:,i) = repmat(gbest(:,:,i),G_p_n,1);      % 扩展leader行向量成为一个G_p_n行的矩阵
            % find visit group 
            F_S = F_S + gbesteval(:,:,i);
            F_V(1,i + 1) = F_S;
            
        end

        % check iter num if can visit 
        % choose visit group
        
        if mod(j1,C_Visit) == 0                                                 % 后期将Iter_times改成j1
            handle = F_V(1,G_n+1)*rand;
            v_t = size(find(handle > F_V(:)),1);                                % 寻找会场组，到这都不用改  
            Vg = visiting(g_pbest,gbestid,G_p_n,v_t,G_n,W1,W2,Vg,Vgxy);
        end
        % 检查循环次数，假如到了重构步，则解散最差的组，并且重新从其他组中选个体组成新的小组
        if mod(j1,C_Recon) == 0
            % 首先全部遣散最差组,在随机指定一个新的leader，最后重构这个组
            [g_pbest, g_pbesteval, g_pvel] = reconstruction(gbesteval, g_pbest, g_pbesteval, g_pvel, gbestid, D, G_n, G_p_n, balance1);
       
            [gbesteval, gbestid] = min(g_pbesteval);
        end
    end
   

  
    Min_global = min(min(g_pbesteval));
    Max_global = max(max(g_pbesteval));
    
    if j1 == 3
        gbesteval_b = squeeze(gbesteval);
        [tongji, g_index] = min(gbesteval_b);
        path = sprintf('%s%d%s%d.mat','FACIPSO_1_3_',g_index,"_",tongji);
        save(path,'g_pbest');
    end
    if j1 == 10
        gbesteval_b = squeeze(gbesteval);
        [tongji, g_index] = min(gbesteval_b);
        path = sprintf('%s%d%s%d.mat','FACIPSO_1_10_',g_index,"_",tongji);
        save(path,'g_pbest');
    end
    if j1 == 20
        gbesteval_b = squeeze(gbesteval);
        [tongji, g_index] = min(gbesteval_b);
        path = sprintf('%s%d%s%d.mat','FACIPSO_1_20_',g_index,"_",tongji);
        save(path,'g_pbest');
    end
    if j1 == 30
        gbesteval_b = squeeze(gbesteval);
        [tongji, g_index] = min(gbesteval_b);
        path = sprintf('%s%d%s%d.mat','FACIPSO_1_30_',g_index,"_",tongji);
        save(path,'g_pbest');
    end
    if j1 == 40
        gbesteval_b = squeeze(gbesteval);
        [tongji, g_index] = min(gbesteval_b);
        path = sprintf('%s%d%s%d.mat','FACIPSO_1_40_',g_index,"_",tongji);
        save(path,'g_pbest');
    end
    if j1 == 50
        gbesteval_b = squeeze(gbesteval);
        [tongji, g_index] = min(gbesteval_b);
        path = sprintf('%s%d%s%d.mat','FACIPSO_1_50_',g_index,"_",tongji);
        save(path,'g_pbest');
    end
    if j1 == 60
        gbesteval_b = squeeze(gbesteval);
        [tongji, g_index] = min(gbesteval_b);
        path = sprintf('%s%d%s%d.mat','FACIPSO_1_60_',g_index,"_",tongji);
        save(path,'g_pbest');
    end
    if j1 == 70
        gbesteval_b = squeeze(gbesteval);
        [tongji, g_index] = min(gbesteval_b);
        path = sprintf('%s%d%s%d.mat','FACIPSO_1_70_',g_index,"_",tongji);
        save(path,'g_pbest');
    end
    if j1 == 80
        gbesteval_b = squeeze(gbesteval);
        [tongji, g_index] = min(gbesteval_b);
        path = sprintf('%s%d%s%d.mat','FACIPSO_1_80_',g_index,"_",tongji);
        save(path,'g_pbest');
    end
    if j1 == 90
        gbesteval_b = squeeze(gbesteval);
        [tongji, g_index] = min(gbesteval_b);
        path = sprintf('%s%d%s%d.mat','FACIPSO_1_90_',g_index,"_",tongji);
        save(path,'g_pbest');
    end
    if j1 == 100
        gbesteval_b = squeeze(gbesteval);
        [tongji, g_index] = min(gbesteval_b);
        path = sprintf('%s%d%s%d.mat','FACIPSO_1_100_',g_index,"_",tongji);
        save(path,'g_pbest');
    end
    if j1 == 150
        gbesteval_b = squeeze(gbesteval);
        [tongji, g_index] = min(gbesteval_b);
        path = sprintf('%s%d%s%d.mat','FACIPSO_1_150_',g_index,"_",tongji);
        save(path,'g_pbest');
    end
    if j1 == 200
        gbesteval_b = squeeze(gbesteval);
        [tongji, g_index] = min(gbesteval_b);
        path = sprintf('%s%d%s%d.mat','FACIPSO_1_200_',g_index,"_",tongji);
        save(path,'g_pbest');
    end
    if j1 == 300
        gbesteval_b = squeeze(gbesteval);
        [tongji, g_index] = min(gbesteval_b);
        path = sprintf('%s%d%s%d.mat','FACIPSO_1_300_',g_index,"_",tongji);
        save(path,'g_pbest');
    end
    if j1 == 400
        gbesteval_b = squeeze(gbesteval);
        [tongji, g_index] = min(gbesteval_b);
        path = sprintf('%s%d%s%d.mat','FACIPSO_1_400_',g_index,"_",tongji);
        save(path,'g_pbest');
    end
    if j1 == 500
        gbesteval_b = squeeze(gbesteval);
        [tongji, g_index] = min(gbesteval_b);
        path = sprintf('%s%d%s%d.mat','FACIPSO_1_500_',g_index,"_",tongji);
        save(path,'g_pbest');
    end
    if j1 == 600
        gbesteval_b = squeeze(gbesteval);
        [tongji, g_index] = min(gbesteval_b);
        path = sprintf('%s%d%s%d.mat','FACIPSO_1_600_',g_index,"_",tongji);
        save(path,'g_pbest');
    end
    if j1 == 700
        gbesteval_b = squeeze(gbesteval);
        [tongji, g_index] = min(gbesteval_b);
        path = sprintf('%s%d%s%d.mat','FACIPSO_1_700_',g_index,"_",tongji);
        save(path,'g_pbest');
    end
    if j1 == 800
        gbesteval_b = squeeze(gbesteval);
        [tongji, g_index] = min(gbesteval_b);
        path = sprintf('%s%d%s%d.mat','FACIPSO_1_800_',g_index,"_",tongji);
        save(path,'g_pbest');
    end
    if j1 == 900
        gbesteval_b = squeeze(gbesteval);
        [tongji, g_index] = min(gbesteval_b);
        path = sprintf('%s%d%s%d.mat','FACIPSO_1_900_',g_index,"_",tongji);
        save(path,'g_pbest');
    end
    if j1 == 1000
        gbesteval_b = squeeze(gbesteval);
        [tongji, g_index] = min(gbesteval_b);
        path = sprintf('%s%d%s%d.mat','FACIPSO_1_1000_',g_index,"_",tongji);
        save(path,'g_pbest');
    end
    if j1 == 1500
        gbesteval_b = squeeze(gbesteval);
        [tongji, g_index] = min(gbesteval_b);
        path = sprintf('%s%d%s%d.mat','FACIPSO_1_1500_',g_index,"_",tongji);
        save(path,'g_pbest');
    end

end
Iter = j1;
gbesteval_b = squeeze(gbesteval);
ii = 1:G_n;
[gbesteval, gbestid] = min(g_pbesteval);
% gbest_b = g_pbest(gbestid(:,:,ii),:,ii);
gbest_b = g_pbest;
[tongji, g_index] = min(gbesteval_b);


%% 头领访问过程(visit processing)
% 学习向量的
% 下面将要把所有的访问过程更新的操作转变成为向量，每个组都使用会场组中的最佳k个个体，暂时选前25%
function visit = visiting(g_pbest,gbestid,G_p_n,v_t,G_n,W1,W2,Vg,Vgxy)

    for i = 1:G_n
        %%%%%%  bsxfun（A,B）的目的是让两个矩阵A，B之间的操作，@times是矩阵相乘，@plus是矩阵相加 @minus是矩阵相减
        Vgxy(:,:,i) = bsxfun(@plus, bsxfun(@times, W1, Vgxy(:,:,i)), bsxfun(@times, W2, bsxfun(@minus, g_pbest(gbestid(:,:,v_t),:,v_t), g_pbest(gbestid(:,:,i),:,i))));  
    end
    Vg = Vgxy;% remove scale parameters
    visit = Vg;
end


%% 小组重构过程(reconstruction processing)
function [g_pbest, g_pbesteval, g_pvel] = reconstruction(gbesteval, g_pbest, g_pbesteval, g_pvel, gbestid, D, G_n, G_p_n, balance1)
    
    [~,I] = min(squeeze(gbesteval));
    
    [g_pbest, g_pbesteval, g_pvel] = releasing(g_pbest, g_pbesteval, g_pvel, D, I, G_n);
    
    % 接下来就是从整体的个体组中随机择新组的leader
    [g_pbest, g_pbesteval, g_pvel] = choosing(g_pbest, g_pbesteval, g_pvel, gbestid, D, I, G_n, balance1);
    % 在接下来是从其他的组里分配给新建组
    [g_pbest, g_pbesteval, g_pvel] = assigning(g_pbest, g_pbesteval, g_pvel, gbestid, I, G_n, G_p_n, balance1);
end


%% tool function area

function [g_pbest, g_pbesteval, g_pvel,fitcount] = DE_M_C_S(g_pbest,g_pbesteval,g_pvel,G_n,F,VRmax,VRmin,func_num,fitcount)
    for jx = 1:G_n
        g_p = g_pbest(:,:,jx);          g_p = rmmissing(g_p);

        tmp_m = g_pbest(:,:,jx);        tmp_m = rmmissing(tmp_m);
        tmp2 = g_pbesteval(:,:,jx);     tmp2 = rmmissing(tmp2);

        g_v = g_pvel(:,:,jx);          g_v = rmmissing(g_v);
        tmp3 = g_pvel(:,:,jx);          tmp3 = rmmissing(tmp3);
        tmp4 = 0.001 .* ones(size(g_pvel,1), size(g_pvel, 2));
        
        NP = size(g_p ,1);
        tab = 1:NP; tab = tab(ones(1,NP),:)';
        dig = 1:NP; DI = (dig - 1) * NP + (1:NP);
        tab(DI) = [];
        tab = reshape(tab, NP-1, [])';
        rand1 = floor(rand(NP,1) * (NP - 1)) + 1;
        rand2 = floor(rand(NP,1) * (NP - 2)) + 2;
        rand3 = floor(rand(NP,1) * (NP - 3)) + 3;
        RND1 = (rand1 - 1) * NP + (1:NP)';
        RND2 = (rand2 - 1) * NP + (1:NP)';
        RND3 = (rand3 - 1) * NP + (1:NP)';
        r1 = tab(RND1); tab(RND1) = tab(:,1);
        % 这两个生成差分向量的从其他小组中随机抽取
        r2 = tab(RND2); tab(RND2) = tab(:,2);
        r3 = tab(RND3);
        % R2 = ceil((G_n - 1) * rand(size(RND1,1),2));
        % R2 = (R2 >= jx) .* (R2 + 1) + (R2 < jx) .* R2;                                    % 抽取小组选取
        % for i = 1:size(RND1,1)
        %     A = g_pbest(:,:,R2(i,1));       B = g_pbest(:,:,R2(i,2));
        %     R3(i,:) = [ceil(size(find(isfinite(sum(A, 2)) & any(A ~= 0, 2)), 1) * rand)...
        %         , ceil(size(find(isfinite(sum(B, 2)) & any(B ~= 0, 2)), 1) * rand)];
        % end
        % r2 = R3(:,1);
        % r3 = R3(:,2);
        % f1 f2 mutation g_pbest(R3(:,1),:,R2(:,1))
        dv = g_p(r1,:) + F .* (g_p(r2,:) - g_p(r3,:));

        dv(bsxfun(@gt,dv, VRmax)) = VRmax;
        dv(bsxfun(@lt,dv, VRmin)) = VRmin;
        % crossover
        cr_lo = rand(size(g_p)) < CR;
        g_p(cr_lo) = dv(cr_lo);
        g_v(cr_lo) = tmp4(cr_lo);
        % selection
        tmp1 = feval(fhd, g_p', func_num); tmp1 = transpose(tmp1); 
        tmp_log = tmp2 > tmp1 & tmp1 > 1; 
        tmp2(tmp_log) = tmp1(tmp_log);
        tmp_m(tmp_log,:) = g_p(tmp_log,:);
        tmp3(tmp_log,:) = g_v(tmp_log,:);
        fitcount = fitcount + size(g_p, 1);
        g_pbest(:,:,jx) = [tmp_m; NaN(size(g_pbest,1) - NP, size(g_pbest, 2))];
        g_pbesteval(:,:,jx) = [tmp2; NaN(size(g_pbesteval,1) - size(tmp2,1),1)];
        g_pvel(:,:,jx) = [tmp3; NaN(size(g_pvel,1) - size(g_v,1), size(g_pvel,2))];
        
    end
    

end

function [g_pbest,g_pbesteval,g_pvel,fitcount] = select(g_pbest_tmp,g_pvel_tmp,g_pbest,g_pbesteval,g_pvel,G_n,func_num,fitcount)
%myFun - Description
%
% Syntax: [g_pbest,g_pbesteval,g_pvel] = select(g_pbest,g_pbesteval,g_pvel,G_n,func_num,fitcount)
%
% Long description
    for jx = 1:G_n
        g_p = g_pbest(:,:,jx);          g_p = rmmissing(g_p);

        tmp_m = g_pbest_tmp(:,:,jx);    tmp_m = rmmissing(tmp_m);
        tmp2 = g_pbesteval(:,:,jx);     tmp2 = rmmissing(tmp2);

        g_v = g_pvel(:,:,jx);          g_v = rmmissing(g_v);
        tmp3 = g_pvel_tmp(:,:,jx);     tmp3 = rmmissing(tmp3);
        tmp1 = feval(fhd, g_p', func_num); tmp1 = transpose(tmp1); 
        tmp_log = tmp2 > tmp1 & tmp1 > 1; 
        tmp2(tmp_log) = tmp1(tmp_log);
        tmp_m(tmp_log,:) = g_p(tmp_log,:);
        tmp3(tmp_log,:) = g_v(tmp_log,:);
        fitcount = fitcount + size(g_p, 1);
        g_pbest(:,:,jx) = [tmp_m; NaN(size(g_pbest,1) - size(tmp_m,1), size(g_pbest, 2))];
        g_pbesteval(:,:,jx) = [tmp2; NaN(size(g_pbesteval,1) - size(tmp2,1),1)];
        g_pvel(:,:,jx) = [tmp3; NaN(size(g_pvel,1) - size(g_v,1), size(g_pvel,2))];
        
    end
end

function [g_pbest, g_pbesteval, g_pvel] = global_restart(VRmin, VRmax, Pop_size, D, Vmin, Vmax, G_n, G_p_n)
%myFun - Description  global_restart function
%
% Syntax: [g_pbest, g_pbesteval, g_pvel] = global_restart(g_pbest, g_pbesteval, g_pvel)
%  global_restart function
% Long description
    pos = bsxfun(@times, VRmin + (VRmax - VRmin), rand(Pop_size, D));
    e = ones(1, Pop_size);
    vel = Vmin + 2 .* Vmax .* rand(Pop_size, D);

    modd = mod(Pop_size, G_n);
    if modd ~= 0
        Pop_size = Pop_size + G_n - modd;
        pos = [pos; NaN((G_n - modd), D)];
        e = [e, NaN(1,(G_n-modd))];
        vel = [vel;NaN((G_n-modd),D)];
    end

    % start page team
    g_pbest = permute(reshape(pos', D, G_p_n, G_n), [2 1 3]);
    g_pvel = permute(reshape(vel', D, G_p_n, G_n), [2 1 3]);
    g_pbesteval = reshape(e', G_p_n, 1, G_n);
end
% local restart
function [g_pbest, g_pbesteval, g_pvel] = local_restart(g_pbest, g_pbesteval, g_pvel)
%myFun - local restart function
%
% Syntax: [g_pbest, g_pbesteval, g_pvel] = local_restart(g_pbest, g_pbesteval, g_pvel)
% 
% Long description
    
end

% choose new leader from other group
function [g_pbest, g_pbesteval, g_pvel] = choosing(g_pbest, g_pbesteval, g_pvel, gbestid, D, I, G_n,balance1)
    % 
    while true
    
        tmp = ceil((G_n - 1)*rand);
        tmp(tmp>=I) = tmp(tmp>=I) + 1;
        if size(find(isfinite(sum(g_pbest(:,:,tmp),2))),1) > balance1
            break;
        end
    end
    A = g_pbest(:,:,tmp);
    tmp1 = ceil((size(find(isfinite(sum(g_pbest(:,:,tmp),2)) & any(A ~= 0, 2)),1)-1)*rand);   
    tmp1= (tmp1 >= gbestid(:,:,tmp)) * (tmp1+1) + (tmp1 < gbestid(:,:,tmp)) * tmp1;

    g_pbest(1,:,I) = g_pbest(tmp1,:,tmp);
    
    A = g_pbest(:,:,tmp);
    A(tmp1,:) = [];
    A = [A;NaN(1,D)];
    g_pbest(:,:,tmp) = A;

    g_pbesteval(1,:,I) = g_pbesteval(tmp1,:,tmp);
    
    A = g_pbesteval(:,:,tmp);
    A(tmp1,:) = [];
    A = [A;NaN(1,1)];
    g_pbesteval(:,:,tmp) = A;           

    g_pvel(1,:,I) = g_pvel(tmp1,:,tmp);
    
    A = g_pvel(:,:,tmp);
    A(tmp1,:) = [];
    A = [A;NaN(1,D)];
    g_pvel(:,:,tmp) = A;
    
end



function [g_pbest, g_pbesteval, g_pvel] = releasing(g_pbest, g_pbesteval, g_pvel, D,  I, G_n)

    A = g_pbest(:,:,I);
    v = size(find(isfinite(sum(g_pbest(:,:,I),2)) & any(A ~= 0, 2)),1);
    % v = size(find(sum(g_pbest(:,:,I),2)~=0),1);         
    tmp = ceil((size(g_pbest, 3)-1)*rand(1,v));                      % 选择解散个体转移的目标组
    tmp(tmp >= I) = tmp(tmp >= I) + 1;
    
    for i = 1:v
        A = g_pbest(:,:,tmp(i));
        g_pbest(size(find(isfinite(sum(g_pbest(:,:,tmp(i)),2)) & any(A ~= 0, 2)),1)+1,:,tmp(i)) = g_pbest(i,:,I);
        % g_pbest(g_pbest == 0) = NaN;
        A = g_pbesteval(:,:,tmp(i));
        g_pbesteval(size(find(isfinite(sum(g_pbesteval(:,:,tmp(i)),2)) & any(A ~= 0, 2)),1)+1,:,tmp(i)) = g_pbesteval(i,:,I);
        % g_pbesteval(g_pbesteval == 0) = NaN;
        A = g_pvel(:,:,tmp(i));
        g_pvel(size(find(isfinite(sum(g_pvel(:,:,tmp(i)),2)) & any(A ~= 0, 2)),1)+1,:,tmp(i)) = g_pvel(i,:,I);
        % g_pvel(g_pvel == 0) = NaN;
        
        g_pbest(i,:,I) = NaN(1,D);
        g_pbesteval(i,:,I) = NaN(1,1);
        g_pvel(i,:,I) = NaN(1,D);
    end
    [g_pbest,g_pbesteval,g_pvel] = myNaN(g_pbest,g_pbesteval,g_pvel,G_n);
end

% 分配成员到新建组
function [g_pbest, g_pbesteval, g_pvel] = assigning(g_pbest, g_pbesteval, g_pvel, gbestid, I, G_n, G_p_n, balance1)
%     k2 = G_p_n-1;
 
%     tmp = ceil((G_n-1)*rand(1,k2));
%     tmp = (tmp >= I) .* bsxfun(@plus, tmp, 1) + (tmp < I) .* tmp;
    % tmp(tmp>=I) = tmp(tmp>=I) + 1;
    i = 1;
    while i < G_p_n
        tmp = ceil((G_n-1)*rand);
        tmp = (tmp >= I) .* (tmp + 1) + (tmp < I) .* tmp;
        if size(find(isfinite(sum(g_pbest(:,:,tmp),2))),1) < balance1
            continue;
        end
        % 从这些组里面提取个体，避免�?�中leader
        A = g_pbest(:,:,tmp);
        % tmp1 = ceil((size(find(isfinite(sum(g_pbest(:,:,tmp),2)) & sum(g_pbest(:,:,tmp),2) ~= 0),1)-1)*rand);
        tmp1 = ceil((size(find(isfinite(sum(g_pbest(:,:,tmp),2)) & any(A ~= 0, 2)),1)-1)*rand);
        if tmp1 == 0
            continue;
        end
        tmp1= (tmp1 >= gbestid(:,:,tmp)) * (tmp1+1) + (tmp1 < gbestid(:,:,tmp)) * tmp1;

        g_pbest(i+1,:,I) = g_pbest(tmp1,:,tmp);
        % g_pbest(g_pbest == 0) = NaN;
        A = g_pbest(:,:,tmp);
        A(tmp1,:) = [];
        A = [A;NaN(1,D)];
        g_pbest(:,:,tmp) = A;

        g_pbesteval(i+1,:,I) = g_pbesteval(tmp1,:,tmp);
        % g_pbesteval(g_pbesteval == 0) = NaN;
        A = g_pbesteval(:,:,tmp);
        A(tmp1,:) = [];
        A = [A;NaN(1,1)];
        g_pbesteval(:,:,tmp) = A;           % 适应度矩阵应该是个列向量

        g_pvel(i+1,:,I) = g_pvel(tmp1,:,tmp);
        % g_pvel(g_pvel == 0) = NaN;
        A = g_pvel(:,:,tmp);
        A(tmp1,:) = [];
        A = [A;NaN(1,D)];
        g_pvel(:,:,tmp) = A;
        i = i + 1;
        
    end
    
end

%% 跳槽过程(hopping processing)
function [pos, e, vel] = hopping(D, pos, e, vel, gbestid, G_n, C_hop, balance1)
    
    C_hop = repmat(C_hop,1,G_n);
    tmp1 = bsxfun(@lt,rand(1,G_n),C_hop);       % 逻辑数组，判断有哪些组的组员可以跳槽
    tmp = ceil((G_n-1) * rand(1,G_n));          % 目标数组粗筛
    tmp2 = 1:G_n;                               % 跳槽源组粗筛
    tmp(tmp1 == 0) = NaN;                       % 允许跳槽的组，�?�择跳槽到的目标
    tmp2(isnan(tmp)) = NaN;                     % 根据目标组，回头确定哪些组是可以跳槽
    tmp2 = rmmissing(tmp2,2);
    tmp = rmmissing(tmp,2);
    tmp3 = zeros(1,size(tmp2,2));               
    
    tmp(tmp >= tmp2) = tmp(tmp >= tmp2) + 1;    % 跳槽目标组的后续步骤&跳槽目标精筛
    
    for i = 1:size(tmp2,2)
        if size(find(isfinite(sum(pos(:,:,tmp2(i)),2))),1) < balance1
            continue;
        end
        
        tmp5 = gbestid(1,1,tmp2(1,i));
        
        tmp6 = rmmissing(pos(:,:,tmp2(i)));
        tmp4 = ceil((size(tmp6,1) -1) * rand);                                  % choose random 
        % tmp4 = ceil((size(find(isfinite(sum(pos(:,:,tmp2(i)),2))),1)-1)*rand);             % tmp the goal is for generate a random number
        if tmp4 == 0
            continue;
        end
        tmp3(1,i) = (tmp4 < tmp5) * tmp4 + (tmp4 >= tmp5) * (tmp4 + 1);             
%         sum(pos(:,:,tmp(i)),2)
%         find(sum(pos(:,:,tmp(i)),2)~=0)

        A = pos(:,:,tmp(i));
        % pos(size(find(isfinite(sum(pos(:,:,tmp(i)),2)) & sum(pos(:,:,tmp(i)),2) ~= 0),1)+1, :, tmp(i)) = pos(tmp3(i),:,tmp2(i));         % success update pos matrix
        pos(size(find(isfinite(sum(pos(:,:,tmp(i)),2)) & any(A ~= 0, 2)),1)+1, :, tmp(i)) = pos(tmp3(i),:,tmp2(i));         % success update pos matrix
        % pos(pos == 0) = NaN;
        B = pos(:,:,tmp2(i));
        B(tmp3(i),:) = [];
        B = [B;NaN(1,D)];
        pos(:,:,tmp2(i)) = B;
        
        A = e(:,:,tmp(i));
        e(size(find(isfinite(sum(e(:,:,tmp(i)),2)) & any(A ~= 0, 2)),1)+1, :, tmp(i)) = e(tmp3(i),:,tmp2(i));             % success update fit matrix
        % e(e == 0) = NaN;
        B = e(:,:,tmp2(i));
        B(tmp3(i),:) = [];
        B = [B;NaN(1,1)];                                                         % add a zero num to the fit matrix
        e(:,:,tmp2(i)) = B;

        A = vel(:,:,tmp(i));
        vel(size(find(isfinite(sum(vel(:,:,tmp(i)),2)) & any(A ~= 0, 2)),1)+1, :, tmp(i)) = vel(tmp3(i),:,tmp2(i));         % success update vel matrix
        % vel(vel == 0) = NaN;
        B = vel(:,:,tmp2(i));
        B(tmp3(i),:) = [];
        B = [B;NaN(1,D)];
        vel(:,:,tmp2(i)) = B;
    end
    % pos(isnan(pos)) = 0;
    % e(isnan(e)) = 0;
    % vel(isnan(vel)) = 0;
    [pos,e,vel] = myNaN(pos,e,vel,G_n);
end
function [g_pbest,g_pbesteval,g_pvel] = myNaN(g_pbest,g_pbesteval,g_pvel,G_n)
%myFun - Description
%
% Syntax: [g_pbest,g_pbesteval,g_pvel] = myNaN(g_pbest,g_pbesteval,g_pvel,G_n)
%
% Long description
    for j = 1:G_n
        A = g_pbest(:,:,j);
        g_pbest(find(all(A == 0,2)),:,j) = NaN(size(find(all(A == 0,2)),1), size(g_pbest,2));
        B = g_pbesteval(:,:,j);
        g_pbesteval(find(all(B == 0,2)),:,j) = NaN(size(find(all(B == 0,2)),1), size(g_pbesteval,2));
        C = g_pvel(:,:,j);
        g_pvel(find(all(C == 0,2)),:,j) = NaN(size(find(all(C == 0,2)),1), size(g_pvel,2));
    end
end
function count = countool(g_pbest, g_pbesteval, g_pvel, G_n)
%myFun - Description
%
% Syntax: countool = countool(g_pbest, g_pbesteval, g_pvel)
%
% Long description
count1 = 0; count2 = 0; count3 = 0; count = 0;
    for i = 1:G_n
        count1 = count1 + size(rmmissing(g_pbest(:,:,i)),1);
        count2 = count2 + size(rmmissing(g_pbesteval(:,:,i)),1);
        count3 = count3 + size(rmmissing(g_pvel(:,:,i)),1);
    end
    count1 = floor(count1/1000);
    count2 = floor(count2/1000);
    count3 = floor(count3/1000);
    count = count1 * 2^2 + count2 * 2^1 + count3 * 2^0;
end
end
