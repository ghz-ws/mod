clear
clc
clf

Fs=10000;   //Sampling frequency
L=10000;    //Length of signal
fc=1000;    //carrier freq
Baud=200;   //data rate
PN=7

//gen. baseband
t=0:1/Fs:(L-1)/Fs;  //time vector. 0-1
lfsr_ask=ones(1,PN);
lfsr_qpsk=ones(1,PN);
lfsr_qam=ones(1,PN);
prbs_tap=[2,3,4,4,6,7,1,6,8,10,1,1,1,15,1,15,12,1,18];  //tap position
rate=Fs/Baud;   //samples per symbol
for i=1:L/rate
    seed_ask=bitxor(lfsr_ask(PN),lfsr_ask(prbs_tap(PN-1)-1));
    lfsr_ask=cat(2,seed_ask,lfsr_ask(1:PN-1));
    ask_dat=seed_ask;
    for k=1:2
        seed_qpsk=bitxor(lfsr_qpsk(PN),lfsr_qpsk(prbs_tap(PN-1)-1));
        lfsr_qpsk=cat(2,seed_qpsk,lfsr_qpsk(1:PN-1));
        qpsk_dat(k)=seed_qpsk;
    end
    for k=1:4
        seed_qam=bitxor(lfsr_qam(PN),lfsr_qam(prbs_tap(PN-1)-1));
        lfsr_qam=cat(2,seed_qam,lfsr_qam(1:PN-1));
        dat_qam(k)=seed_qam;
    end
    for j=1:rate
        ask_bb(rate*(i-1)+j)=ask_dat;             //format 1/0
        bpsk_bb(rate*(i-1)+j)=(ask_dat-0.5)*2;    //format 1/-1
        qpsk_ibb(rate*(i-1)+j)=(qpsk_dat(1)-0.5)*2;   //format 1/-1
        qpsk_qbb(rate*(i-1)+j)=(qpsk_dat(2)-0.5)*2;   //format 1/-1
        qam_ibb(rate*(i-1)+j)=(dat_qam(1)+2*dat_qam(2)-1.5)/1.5;    //format 1/0.3/-0.3/-1
        qam_qbb(rate*(i-1)+j)=(dat_qam(3)+2*dat_qam(4)-1.5)/1.5;    //format 1/0.3/-0.3/-1
    end
end

//modulation
wave1=cos(2*%pi*fc*t);
wave2=sin(2*%pi*fc*t);
for i=1:L
    ask(1,i)=ask_bb(i)*wave1(i);
    bpsk(1,i)=bpsk_bb(i)*wave1(i);
    qpsk(1,i)=qpsk_ibb(i)*wave1(i)+qpsk_qbb(i)*wave2(i);
    qam(1,i)=qam_ibb(i)*wave1(i)+qam_qbb(i)*wave2(i);
end

//fft
f=Fs*(0:(L/2))/L; //freq vector. freq res.=Fs/N*n. n=0-N/2
fn=length(f) //length of freq
win=window('hn',L)      //window func.
y1(1,:)=fft(ask.*win);
y2(1,:)=fft(bpsk.*win);
y3(1,:)=fft(qpsk.*win);
y4(1,:)=fft(qam.*win);

subplot(4,2,1)
plot(t,ask,"thickness",2)
xgrid()
title('ASK Waveform',"fontsize",5)
xlabel('Time [s]',"fontsize",5)
ylabel('Amplitude',"fontsize",5)
g1=gca();   //get axis object
g1.data_bounds(:,1)=[0;1/Baud*10];    //x axis scale
g1.font_size=5

subplot(4,2,3)
plot(t,bpsk,"thickness",2)
xgrid()
title('BPSK Waveform',"fontsize",5)
xlabel('Time [s]',"fontsize",5)
ylabel('Amplitude',"fontsize",5)
g2=gca();   //get axis object
g2.data_bounds(:,1)=[0;1/Baud*10];    //x axis scale
g2.font_size=5

subplot(4,2,5)
plot(t,qpsk,"thickness",2)
xgrid()
title('QPSK Waveform',"fontsize",5)
xlabel('Time [s]',"fontsize",5)
ylabel('Amplitude',"fontsize",5)
g3=gca();   //get axis object
g3.data_bounds(:,1)=[0;1/Baud*10];    //x axis scale
g3.font_size=5

subplot(4,2,7)
plot(t,qam,"thickness",2)
xgrid()
title('16QAM Waveform',"fontsize",5)
xlabel('Time [s]',"fontsize",5)
ylabel('Amplitude',"fontsize",5)
g4=gca();   //get axis object
g4.data_bounds(:,1)=[0;1/Baud*10];    //x axis scale
g4.font_size=5

//show spectrum
subplot(4,2,2)
plot(f,20*log10(abs(y1(1:fn))),"thickness",2)
xgrid()
title("ASK Spectrum","fontsize",5)
xlabel("Freq [Hz]","fontsize",5)
ylabel("Amplitude [dB]","fontsize",5)
g5=gca();   //get axis object
g5.data_bounds(:,1)=[0;fc*2];    //x axis scale
g5.data_bounds(:,2)=[0;60];    //y axis scale
g5.font_size=5

subplot(4,2,4)
plot(f,20*log10((abs(y2(1:fn)))),"thickness",2)
xgrid()
title("BPSK Spectrum","fontsize",5)
xlabel("Freq [Hz]","fontsize",5)
ylabel("Amplitude [dB]","fontsize",5)
g6=gca();   //get axis object
g6.data_bounds(:,1)=[0;fc*2];    //x axis scale
g6.data_bounds(:,2)=[0;60];    //y axis scale
g6.font_size=5

subplot(4,2,6)
plot(f,20*log10((abs(y3(1:fn)))),"thickness",2)
xgrid()
title("QPSK Spectrum","fontsize",5)
xlabel("Freq [Hz]","fontsize",5)
ylabel("Amplitude [dB]","fontsize",5)
g7=gca();   //get axis object
g7.data_bounds(:,1)=[0;fc*2];    //x axis scale
g7.data_bounds(:,2)=[0;60];    //y axis scale
g7.font_size=5

subplot(4,2,8)
plot(f,20*log10((abs(y4(1:fn)))),"thickness",2)
xgrid()
title("16QAM Spectrum","fontsize",5)
xlabel("Freq [Hz]","fontsize",5)
ylabel("Amplitude [dB]","fontsize",5)
g8=gca();   //get axis object
g8.data_bounds(:,1)=[0;fc*2];    //x axis scale
g8.data_bounds(:,2)=[0;60];    //y axis scale
g8.font_size=5
