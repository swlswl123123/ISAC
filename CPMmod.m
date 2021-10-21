function [CPM_BB] = CPMmod(data, oversamp)
    %myFun - Description
    %
    % Syntax: [output] = myFun(input)
    %
    % Long description
    
    bit_rate = 10;    % 符号速率MHz
    Tb = 1/bit_rate;  % 符号时间us   
    fs = bit_rate * oversamp; % 采样率
    dt = 1/fs;        % 采样间隔us
    ll = 2;           % 关联长度
    
    %% 生成高斯滤波函数gaussf与其积分表示
    tt = -Tb : 1/fs: Tb - 1/fs;
    tt = tt + 1/fs/2;
    gaussf = ones(size(tt));
    % plot(gaussf)
    for i = 1:length(gaussf)
    
        if i == 1
            J_g(i) = 0;
        else
            J_g(i) = J_g(i - 1) + gaussf(i) * dt;
        end
    end
    J_g = J_g / 4 / Tb;
    % figure
    % plot(J_g)
    
    %% 调制MSK
    data = 2*data - 1;
    % for i = 1:length(data)
    %     if i > 1
    %         data(i) = data(i) * data(i-1);
    %     end
    % end
    bit = [];
    phi_all = [];
    L = 0;
    phi = zeros(1, oversamp);
    for i = 1:length(data)
        if i == 1
            if data(1) == 1
                phi_all = [phi_all, ((1:oversamp)-1)*pi/2/oversamp];
            else
                phi_all = [phi_all, -((1:oversamp)-1)*pi/2/oversamp];
            end
        else
            bit = data(i-1:i);
            if i > 2
                L = L + data(i-2);
            end
            for j  = 1:oversamp
                phi(j) = bit(2)*J_g(j) + bit(1)*J_g(j + oversamp);
            end
            phi_all = [phi_all, pi * phi + L * pi / 2 + data(1) * pi / 4];
        end
    end
    phi_all = mod(phi_all, 2*pi);
    % figure
    % t = 0:dt:length(data)*Tb-dt;
    % plot(mod(phi_all, 2*pi))
    CPM_BB = exp(1i*phi_all);
    end