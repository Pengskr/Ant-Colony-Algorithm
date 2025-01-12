clc
clear
close all

%% 初始化
% 根据节点的邻近节点表及字母节点-数字节点对应表，构造节点元胞数组
nodes_data = cell(0);
nodes_data(1, :) = {1, [2, 6, 7], [12, 16, 14]};
nodes_data(2, :) = {2, [1, 3, 6], [12, 10, 7]};
nodes_data(3, :) = {3, [2, 4, 5, 6], [10, 3, 5, 6]};
nodes_data(4, :) = {4, [3, 5], [3, 4]};
nodes_data(5, :) = {5, [3, 4, 6, 7], [5, 4, 2, 8]};
nodes_data(6, :) = {6, [1, 2, 3, 5, 7], [16, 7, 6, 2, 9]};
nodes_data(7, :) = {7, [1, 5, 6], [14, 8, 9]};

% 定义起点和终点
node_start = 4;
node_end = 1;

% 蚁群相关定义
m = 50;                     % 蚁群数量
n = size(nodes_data, 1);    % 节点数量
alpha = 1;                  % 信息素启发式因子α取值范围一般为[l， 4]
beta = 3;                   % 期望启发因子β取值范围一般为[3， 5]
                            % 在后续测试中发现，beta=5，算法容易陷入局部最优，因为在地图中F点时，选择A还是B受到启发因子的影响，会偏向于选择B，因此造成局部最优

rho = 0.1;                  % 信息素蒸发率
Q = 1;                      % 相关人员研究结果表明： 总信息量 Q 对 ant-cycle 模型蚁群算法的性能没有明显的影响。 

iter_max = 100;                         % 一般取100~500
Route_best = cell(iter_max,1);          % 各代最佳路径       
Length_best = zeros(iter_max,1);        % 各代最佳路径的长度  
Length_ave = zeros(iter_max,1);         % 各代路径的平均长度

% 将信息素、启发因子一并放入nodes_data中
Delta_Tau_initial = nodes_data(:,1:2);
for i = 1:n
    nodes_data{i,4} = ones(1, length(nodes_data{i,3}));     % 各路径上的信息素初始化为1
    nodes_data{i,5} = 1./nodes_data{i,3};                   % 启发因子：一般为距离的倒数
    Delta_Tau_initial{i,3} = zeros(1, length(nodes_data{i,3}));
end

%% 迭代寻优
for iter = 1 : iter_max
    
    route = cell(0);    % 存放第iter代每只蚂蚁的路线

    for i = 1 : m
        node_step = node_start;     % 当前所在节点
        path = [node_step, ];       % 记录第i只蚂蚁走过的路径
        dist = 0;                   % 已走过的距离

        while ~ismember(node_end, path) % 当路径表里面包含了终节点时，该蚂蚁完成路径寻优，跳出循环
            % 寻找邻近节点
            neighbor = nodes_data{node_step, 2};

            % 删除邻近节点中已经访问过的节点
            idx = [];
            for k = 1 : length(neighbor)
                if ismember(neighbor(k),path)
                    idx(end+1) = k;
                end
            end
            neighbor(idx) = [];

            % 判断是否进入死胡同， 若是，直接返回到起点，重新寻路
            if isempty(neighbor)
                node_step = node_start;
                path = [node_step, ];
                dist = 0;
                continue
            end

            % 计算下一个节点的访问概率
            P = neighbor;
            for k = 1 : length(P)
                idx = find(nodes_data{node_step, 2} == neighbor(k));
                P(2, k) = nodes_data{node_step, 4}(idx)^alpha * ...  %%%%%
                          nodes_data{node_step, 5}(idx)^beta;
            end
            P(2, :) = P(2, :) / sum(P(2, :));

            % 轮盘赌法 选择下一个访问节点
            Pc = cumsum(P(2, :));   % B = cumsum(A) 从 A 中的第一个其大小不等于 1 的数组维度开始返回 A 的累积和。
            Pc = [0, Pc];
            randnum = rand;
            for k = 1 : length(Pc) - 1
                if randnum > Pc(k) && randnum < Pc(k+1)
                    target_node = neighbor(k);
                end
            end
%             % 概率最大 作为下一个访问节点 (经实验对比，这种方式会丢失掉随机性，结果为一条水平直线)
%             [Pmax, idx_Pmax] = max(P(2, :));
%             target_node = neighbor(idx_Pmax);


            % 计算单步距离
            idx = find(nodes_data{node_step, 2} == target_node);
            dist = dist + nodes_data{node_step, 3}(idx);

            % 更新下一步的目标节点及路径集合
            node_step = target_node;
            path(end + 1) = node_step;
        end
        
        % 存放第i只蚂蚁的累计距离以及对应路径
        Length(i, 1) = dist;
        route{i, 1} = path;
    end

%     % 计算到第iter代为止的最短距离及对应路径
%     if iter == 1
%         [min_Length, min_idx] = min(Length);
%         Length_best(iter) = min_Length;
%         Length_ave(iter) = mean(Length);
%         Route_best{iter, 1} = route{min_idx, 1};
%     else
%         [min_Length, min_idx] = min(Length);
%         Length_best(iter) = min(Length_best(iter-1), min_Length);
%         Length_ave(iter) = mean(Length);
%         if Length_best(iter) == min_Length
%             Route_best{iter, 1} = route{min_idx, 1};
%         else
%             Route_best{iter, 1} = Route_best{iter-1, 1};
%         end
%     end

    [min_Length, min_idx] = min(Length);
    Length_best(iter) = min_Length;
    Length_ave(iter) = mean(Length);
    Route_best{iter, 1} = route{min_idx, 1};


    % 更新信息素
    Delta_Tau = Delta_Tau_initial;
    for i = 1 : m
        for j = 1 : length(route{i, 1})-1
            node_start_temp = route{i,1}(j);
            node_end_temp = route{i,1}(j+1);
            idx =  find(Delta_Tau{node_start_temp, 2} == node_end_temp);
            Delta_Tau{node_start_temp,3}(idx) = Delta_Tau{node_start_temp,3}(idx) + Q/Length(i);
        end
    end

    % 考虑挥发因子，更新信息素
    for i = 1 : n
        nodes_data{i, 4} = (1-rho) * nodes_data{i, 4} + Delta_Tau{i, 3};
    end  
end

%% 绘图、结果    
figure
plot(1 : iter_max, Length_best, 'b', 1 : iter_max, Length_ave, 'r')
legend('最短距离','平均距离')
xlabel('迭代次数')
ylabel('距离')
title('各代最短距离与平均距离对比')

% 最优路径
[dist_min, idx] = min(Length_best);
path_opt = Route_best{idx,1}
