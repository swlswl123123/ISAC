function [CPM_BB] = CPMmod(data, oversamp)
    %myFun - Description
    %
    % Syntax: [output] = myFun(input)
    %
    % Long description
    
    %% 参数
    bit_rate = 10;    % 符号速率MHz
    Tb = 1/bit_rate;  % 符号时间us   
    fs = bit_rate * oversamp; % 采样率
    dt = 1/fs;        % 采样间隔us
    ll = 2;           % 关联长度
    
    %% 调制MSK
    bit = [];
    phi_all = [];
    L = 0;
    phasePre = 0;
    for i = 1:length(data)
        if i == 1
            bit = [0, data(i)];
        else
            bit = data(i-1:i);
        end
        if bit == [1 1]
            phi_all = [phi_all, ((1:oversamp)-1)*pi/2/oversamp + phasePre];
            phasePre = phasePre + pi/2;
        elseif bit == [0 0]
            phi_all = [phi_all, -((1:oversamp)-1)*pi/2/oversamp + phasePre];
            phasePre = phasePre - pi/2;
        else
            phi_all = [phi_all, ones(1, int32(oversamp))*phasePre];
        end
    end
    phi_all = mod(phi_all, 2*pi);
    % figure
    % t = 0:dt:length(data)*Tb-dt;
    % plot(mod(phi_all, 2*pi))
    CPM_BB = exp(1i*phi_all);
    end