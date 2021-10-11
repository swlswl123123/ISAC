function out = MSKmod(data, oversample)
    bitstream = 2*data - 1;
    N=length(bitstream);
    %实现差分编码
    b0=1;
    for i=1:N
        encode_output(i)=b0*bitstream(i);%对应bk
        b0=encode_output(i);
    end
    I=[];Q=[];%奇数进I路、偶数进Q路
    for i=1:N
        if mod(i,2)==0
            I=[I,encode_output(i),encode_output(i)];
        else
            Q=[Q,encode_output(i),encode_output(i)];
        end
    end
    bit_data=[];
    for i=1:N
        bit_data=[bit_data,bitstream(i)*ones(1,oversample)];%脉冲成型，认为输入的比特流是矩形的
    end
    I=[1,I];%把最开始的I1补上
    I=I(1:N);
    Q=Q(1:N);
    Q=-Q;
    I_data=[];
    Q_data=[];
    base_wave = 0 : 1 : oversample - 1;
    for i=1:N
        I_data=[I_data,I(i)*cos(pi*((i-1)+base_wave)/2/oversample)];%I路:
        Q_data=[Q_data,Q(i)*sin(pi*((i-1)+base_wave)/2/oversample)];%Q
    end
    out = I_data + 1i*Q_data;
    % figure();
    % t=0:1/Fs:N*T-1/Fs;
    % subplot(3,1,1)
    % plot(t,bit_data);legend('比特流')
    % subplot(3,1,2)
    % plot(t,I_data);legend('I路比特流')
    % subplot(3,1,3)
    % plot(t,Q_data);legend('Q路比特流')
    %载波信号
    % bit_t=0:1/Fs:N*T-1/Fs;%定义时间轴
    % %定义I路和Q路的载波信号
    % I_carrier=cos(2*pi*(fc+fd)*bit_t);
    % Q_carrier=sin(2*pi*(fc+fd)*bit_t);
    % I_dot=I_data.*I_carrier;
    % Q_dot=Q_data.*Q_carrier;
    % Dot=I_dot+Q_dot;
    % figure();
    % subplot(3,1,1)
    % plot(t,Dot);legend('总信号')
    % subplot(3,1,2)
    % plot(t,I_dot);legend('I路调制后信号')
    % subplot(3,1,3)
    % plot(t,Q_dot);legend('Q路调制后信号')
    % figure();
    % plot(10*log10(abs(fft(Dot))));legend('总信号频谱')
end