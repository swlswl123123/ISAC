function result = GMSK_demod_new(data_store, data_num, inter_wave, mode, u, tag);
%myFun - Description
%
% Syntax: result = GMSK_demod_new(data_store, data_num, inter_wave, mode)
%
% Long description

num_data_frame = 1024;
oversamp_BB_1 = 2;
oversamp_BB_2 = 2;
% load('lib/g_1024.mat');  % GMSK调制 g函数
% g = g(1:16:end);
% cnt = csvread("count.csv");

switch mode
case 1
    num_pulse = 12;
    KK = 3;
case 2
    num_pulse = 12;
    KK = 3;
case 3
    num_pulse = 48;
    KK = 5;
case 4
    num_pulse = 96;
    KK = 5;
end

soft_info_pre = [];
% 10 ~ 15 64阶
LPF = [7.42676310266323e-05,3.85670800492053e-05,-0.000153183665210777,-0.000257210496571457,7.49898328952584e-05,0.000583650771968760,0.000385767086817913,-0.000722526213769458,-0.00131681919121312,7.69135387640725e-05,0.00222723929667286,0.00175894964056264,-0.00201388157764273,-0.00436720864542797,-0.000561095955837434,0.00600868961247023,0.00569016345001960,-0.00408072347733571,-0.0114009302689952,-0.00333995992364471,0.0133837086703058,0.0153518542983937,-0.00651988043147207,-0.0267953601030120,-0.0122935325491800,0.0285059936392965,0.0414523144822695,-0.00853430272531288,-0.0739467168470990,-0.0529700180163418,0.0995770527205206,0.298717740819656,0.390681074632244,0.298717740819656,0.0995770527205206,-0.0529700180163418,-0.0739467168470990,-0.00853430272531288,0.0414523144822695,0.0285059936392965,-0.0122935325491800,-0.0267953601030120,-0.00651988043147207,0.0153518542983937,0.0133837086703058,-0.00333995992364471,-0.0114009302689952,-0.00408072347733571,0.00569016345001960,0.00600868961247023,-0.000561095955837434,-0.00436720864542797,-0.00201388157764273,0.00175894964056264,0.00222723929667286,7.69135387640725e-05,-0.00131681919121312,-0.000722526213769458,0.000385767086817913,0.000583650771968760,7.49898328952584e-05,-0.000257210496571457,-0.000153183665210777,3.85670800492053e-05,7.42676310266323e-05];

if mode == 2
    de_out_all = [];
    for v = 1:num_pulse
        th = (data_num(v) - 256)/2;
        data1 = data_store(v,1:th*oversamp_BB_1);
        % data1 = best_sample_interpreter(data1, u, tag);
        % 插值
        % data1_oversample4 = zeros(1, length(data1)*2);
        % data1_oversample4(1:2:end) = data1;
        % data1_tmp = conv(data1_oversample4, LPF);
        % data1_tmp = data1_tmp(32+1:32+length(data1_oversample4));
        % data1_tmp(1) = abs(data1(1));
        % [~, soft_info] = GMSK_viterbi(data1_tmp, th, oversamp_BB_2, 0, g);
        if mod(th, 2) == 0
            [de_out, soft_info] = GMSK_viterbi_origin(data1, oversamp_BB_2, 0, th);
        else
            [de_out, soft_info] = GMSK_viterbi_origin(data1, oversamp_BB_2, 1, th+1);
        end
        soft_info_pre = [soft_info_pre, soft_info];
        de_out_all = [de_out_all, de_out];

        data2 = data_store(v,1+th*oversamp_BB_1:(th+256)*oversamp_BB_1);
        % data2 = best_sample_interpreter(data2, u, tag);

        % 插值
        % data2_oversample4 = zeros(1, length(data2)*2);
        % data2_oversample4(1:2:end) = data2;
        % data2_tmp = conv(data2_oversample4, LPF); 
        % data2_tmp = data2_tmp(32+1:32+length(data2_oversample4));
        % data2_tmp(1) = abs(data2(1));
        % [~, soft_info] = GMSK_viterbi(data2_tmp, 256, oversamp_BB_2, 0, g);
        [de_out, soft_info] = GMSK_viterbi_origin(data2, oversamp_BB_2, 0, 256);
        soft_info_pre = [soft_info_pre, soft_info];
        de_out_all = [de_out_all, de_out];

        data3 = data_store(v,1+(th+256)*oversamp_BB_1:(2*th+256)*oversamp_BB_1);
        % data3 = best_sample_interpreter(data3, u, tag);
        % if v == 7
        %     plot(angle(data3))
        %     hold off;
        % end
        % 插值
        % data3_oversample4 = zeros(1, length(data3)*2);
        % data3_oversample4(1:2:end) = data3;
        % data3_tmp = conv(data3_oversample4, LPF); 
        % data3_tmp = data3_tmp(32+1:32+length(data3_oversample4));
        % data3_tmp(1) = abs(data3(1));
        % [~, soft_info] = GMSK_viterbi(data3_tmp, th, oversamp_BB_2, 0, g);
        if mod(th, 2) == 0
            [de_out, soft_info] = GMSK_viterbi_origin(data3, oversamp_BB_2, 0, th);
        else
            [de_out, soft_info] = GMSK_viterbi_origin(data3, oversamp_BB_2, 1, th+1);
        end
        soft_info_pre = [soft_info_pre, soft_info];
        de_out_all = [de_out_all, de_out];
    end
    % csvwrite("de_out_all.csv", de_out_all);
else
    de_out_all = [];
    for v = 1:num_pulse
        % viterbi译码
        % 插值
        data = data_store(v,:);
        % data = best_sample_interpreter(data, u, tag);
        % if v == 1
        %     plot(real(data))
        %     hold off;
        % end
        % data_oversample4 = zeros(1, length(data)*2);
        % data_oversample4(1:2:end) = data;
        % data_tmp = conv(data_oversample4, LPF);
        % data_tmp = data_tmp(32+1:32+length(data_oversample4));
        % data_tmp(1) = abs(data(1));
        % [~, soft_info] = GMSK_viterbi(data_tmp, data_num(v), oversamp_BB_2, 0, g);
        [de_out, soft_info] = GMSK_viterbi_origin(data, oversamp_BB_2, 0, data_num(v));
        % 软信息组合
        soft_info_pre = [soft_info_pre, soft_info(36+1:36+256)];
        de_out_all = [de_out_all, de_out];
    end
    % csvwrite("de_out_all" + num2str(cnt) + ".csv", de_out_all);
    % csvwrite("de_out_all.csv", de_out_all);
end

% 解交织
[~, index] = sort(inter_wave);
% soft_info_pre = soft_info_pre(index);

switch mode
case 1
    % csvwrite("soft_info_pre" + num2str(cnt) + ".csv", soft_info_pre);
    soft_info_final = soft_info_pre(index);
case 2
    soft_info_pre_p1 = soft_info_pre(1:1024*3);
    soft_info_pre_p2 = soft_info_pre(1+1024*3:1024*6);
    % csvwrite("soft_info_pre_p1.csv", soft_info_pre_p1);
    % csvwrite("soft_info_pre_p2.csv", soft_info_pre_p2);
    soft_info_final = soft_info_pre_p1(index) + soft_info_pre_p2(index);
case 3
    soft_info_pre_p1 = soft_info_pre(1:1024*5);
    soft_info_pre_p2 = soft_info_pre(1+1024*5:1024*10);
    soft_info_final = soft_info_pre_p1(index) + soft_info_pre_p2(index);
case 4
    soft_info_pre_p1 = soft_info_pre(1:1024*5);
    soft_info_pre_p2 = soft_info_pre(1+1024*5:1024*10);
    soft_info_pre_p3 = soft_info_pre(1+1024*10:1024*15);
    soft_info_pre_p4 = soft_info_pre(1+1024*15:1024*20);
    soft_info_final = soft_info_pre_p1(index) + soft_info_pre_p2(index) + soft_info_pre_p3(index) + soft_info_pre_p4(index);
end

% result = decoderTurbo5(-soft_info_final', interleaver(512));
% cnt = cnt + 1;
% csvwrite("count.csv", cnt);
[result, ~] = Ldecoder2(soft_info_final', num_data_frame, num_data_frame*KK);
result = result';
    
end