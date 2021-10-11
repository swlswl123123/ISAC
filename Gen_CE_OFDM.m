%% 
% 任务要求 100MHz，50us
% 第一次： 
% 设置10个符号，每个符号时长5us，有效时长4us，循环前缀长度1us，共400个有效子载波，512点IFFT，16QAM，信息传输速率400Mbit/s
function frame = Gen_CE_OFDM(Num_symbol_data,Br,T_subframe)
    % Br = 100MHz
    % T_subframe = 5us
    % Num_sysmbol_data = 10
    
    % Br = 1e8;
    % T_subframe = 5e-6;
    % Num_symbol_data = 10;
    B = Br;
    Ts = T_subframe*4/5; %有效时长
    Tg = T_subframe/5;%循环前缀
    T = Ts+Tg;
    carrier_f = 1/Ts;%子载波间隔
    Num_carrier = floor(B*Ts); %有效载波数
    ifft_length = 2^nextpow2(Num_carrier)*2;%IFFT点数
    cp_length = floor(ifft_length/4); 
    Nr = ifft_length + cp_length;
    dt = Ts/ifft_length;
    Fsr = 1/dt;
    modulate_bit = 4;    %1为BPSK，2为QPSK,TD-LTE中用的是16QAM，即为4
    %% 参数输出
    fprintf('带宽：');
    fprintf('Br = %d MHz\n',B/1e6);
    fprintf('采样率：');
    fprintf('fsr = %d MHz\n',int32(Fsr/1e6));
    fprintf('OFDM块时长：');
    fprintf('Tg = %d us\n',int32(T/1e-6));
    fprintf('载波个数：');
    fprintf('Nc = %d\n',Num_carrier);
    %% 生成bit码
    inforSource = randi([0 1],1,Num_carrier*modulate_bit*Num_symbol_data/2);
    % inforSource = ones(1,Num_carrier*modulate_bit*Num_symbol_data);
    % save inforSource inforSource;
    data_temp1 = reshape(inforSource,modulate_bit,[]).';
    lteSymMap = [11 10 14 15 9 8 12 13 1 0 4 5 3 2 6 7];
    %% 星座图映射
    modulate_data = qammod(bi2de(data_temp1), 2^(modulate_bit), lteSymMap, 'UnitAveragePower', true);
    % scatterplot(modulate_data);title('星座图');
    %% 子载波映射
    carrier_data_part1 = reshape(modulate_data,Num_carrier/2,Num_symbol_data);
    idx = fliplr(1:Num_carrier/2);
    carrier_data_part2 = conj(carrier_data_part1(idx, :));
    carrier_data_temp1 = [zeros(1, Num_symbol_data); carrier_data_part1; zeros(ifft_length-Num_carrier-2, Num_symbol_data); zeros(1, Num_symbol_data); carrier_data_part2];
    carrier_data = (carrier_data_temp1);
    % carrier_data = carrier_data_temp1;
    %% OFDM调制
    time_data_ifft = ifft(carrier_data);
    frame = reshape(time_data_ifft,1,[]);
    %% 加CP

    % time_data_cp = [time_data_ifft(end-cp_length+1:end,:);time_data_ifft];
    %% 波形生成
    % frame = reshape(time_data_cp,1,[]);
    %% 画图
    % figure,plot(dt*(0:length(frame)-1),abs(frame));
    % xlabel('时间/s');ylabel('幅度');title('发射波形');
    