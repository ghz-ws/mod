clear
clc
clf

Fs=10000     //sample freq
fc=1000;    //carrier freq
Baud=200    //symbol rate
alpha=0.5   //roll off
TAP=261      //TAP
L=4000;     //Length of signal
PN=7        //PRBS n. n=2-20 except n=8,12,13,14,16,19

//calc rc impulse response
T=1/Baud;
Ts=1/Fs;
for i=1:Fs/Baud*10
    t=i*Ts;
    time(i)=t;
    h(i)=sin(%pi*t/T)*T/(%pi*t)*cos(%pi*t*alpha/T)/(1-(2*t*alpha/T)^2)
end

//extract coefficients
half_tap=int(TAP/2)
coe(half_tap+1)=1;
for i=1:half_tap
    coe(i)=h(half_tap+1-i);
    coe(half_tap+i+1)=h(i);
    if coe(i)==%inf
        coe(i)=0;
    end
    if coe(half_tap+i+1)==%inf
        coe(half_tap+i+1)=0;
    end
end

//gen. time vector and baseband
for i=1:L
    X(i)=(i-1)*1/Fs;    //time vector
end
lfsr_ask=ones(1,PN);
lfsr_bpsk=ones(1,PN);
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
        qam_dat(k)=seed_qam;
    end
    for j=1:rate
        if j==1
            ask_bb(rate*(i-1)+j)=ask_dat;             //format 1/0
            bpsk_bb(rate*(i-1)+j)=(ask_dat-0.5)*2;    //format 1/-1
            qpsk_ibb(rate*(i-1)+j)=(qpsk_dat(1)-0.5)*2;   //format 1/-1
            qpsk_qbb(rate*(i-1)+j)=(qpsk_dat(2)-0.5)*2;   //format 1/-1
            qam_ibb(rate*(i-1)+j)=(qam_dat(1)+2*qam_dat(2)-1.5)/1.5;    //format 1/0.3/-0.3/-1
            qam_qbb(rate*(i-1)+j)=(qam_dat(3)+2*qam_dat(4)-1.5)/1.5;    //format 1/0.3/-0.3/-1
        else
            ask_bb(rate*(i-1)+j)=0;
            bpsk_bb(rate*(i-1)+j)=0;
            qpsk_ibb(rate*(i-1)+j)=0;
            qpsk_qbb(rate*(i-1)+j)=0;
            qam_ibb(rate*(i-1)+j)=0;
            qam_qbb(rate*(i-1)+j)=0;
        end
    end
end

//process filter
lfsr_ask_filter=zeros(1,size(coe)(1));
lfsr_bpsk_filter=zeros(1,size(coe)(1));
lfsr_iqpsk_filter=zeros(1,size(coe)(1));
lfsr_qqpsk_filter=zeros(1,size(coe)(1));
lfsr_iqam_filter=zeros(1,size(coe)(1));
lfsr_qqam_filter=zeros(1,size(coe)(1));
for i=1:L
    lfsr_ask_filter=cat(2,ask_bb(i),lfsr_ask_filter(1:size(coe)(1)-1));
    lfsr_bpsk_filter=cat(2,bpsk_bb(i),lfsr_bpsk_filter(1:size(coe)(1)-1));
    lfsr_iqpsk_filter=cat(2,qpsk_ibb(i),lfsr_iqpsk_filter(1:size(coe)(1)-1));
    lfsr_qqpsk_filter=cat(2,qpsk_qbb(i),lfsr_qqpsk_filter(1:size(coe)(1)-1));
    lfsr_iqam_filter=cat(2,qam_ibb(i),lfsr_iqam_filter(1:size(coe)(1)-1));
    lfsr_qqam_filter=cat(2,qam_qbb(i),lfsr_qqam_filter(1:size(coe)(1)-1));
    ask_bb_filtered(i)=sum(coe'.*lfsr_ask_filter);
    bpsk_bb_filtered(i)=sum(coe'.*lfsr_bpsk_filter);
    qpsk_ibb_filtered(i)=sum(coe'.*lfsr_iqpsk_filter);
    qpsk_qbb_filtered(i)=sum(coe'.*lfsr_qqpsk_filter);
    qam_ibb_filtered(i)=sum(coe'.*lfsr_iqam_filter);
    qam_qbb_filtered(i)=sum(coe'.*lfsr_qqam_filter);
end

tc=0:1/Fs:(L-1)/Fs;  //time vector. 0-1
wave1=cos(2*%pi*fc*tc);
wave2=sin(2*%pi*fc*tc);
for i=1:L
    ask(1,i)=ask_bb_filtered(i)*wave1(i);
    bpsk(1,i)=bpsk_bb_filtered(i)*wave1(i);
    qpsk(1,i)=qpsk_ibb_filtered(i)*wave1(i)+qpsk_qbb_filtered(i)*wave2(i);
    qam(1,i)=qam_ibb_filtered(i)*wave1(i)+qam_qbb_filtered(i)*wave2(i);
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
plot(tc,ask,"thickness",2)
xgrid()
title('ASK Waveform',"fontsize",5)
xlabel('Time [s]',"fontsize",5)
ylabel('Amplitude',"fontsize",5)
g1=gca();   //get axis object
g1.data_bounds(:,1)=[0;1/Baud*10];    //x axis scale
g1.font_size=5

subplot(4,2,3)
plot(tc,bpsk,"thickness",2)
xgrid()
title('BPSK Waveform',"fontsize",5)
xlabel('Time [s]',"fontsize",5)
ylabel('Amplitude',"fontsize",5)
g2=gca();   //get axis object
g2.data_bounds(:,1)=[0;1/Baud*10];    //x axis scale
g2.font_size=5

subplot(4,2,5)
plot(tc,qpsk,"thickness",2)
xgrid()
title('QPSK Waveform',"fontsize",5)
xlabel('Time [s]',"fontsize",5)
ylabel('Amplitude',"fontsize",5)
g3=gca();   //get axis object
g3.data_bounds(:,1)=[0;1/Baud*10];    //x axis scale
g3.font_size=5

subplot(4,2,7)
plot(tc,qam,"thickness",2)
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
g5.data_bounds(:,2)=[-30;50];    //y axis scale
g5.font_size=5

subplot(4,2,4)
plot(f,20*log10((abs(y2(1:fn)))),"thickness",2)
xgrid()
title("BPSK Spectrum","fontsize",5)
xlabel("Freq [Hz]","fontsize",5)
ylabel("Amplitude [dB]","fontsize",5)
g6=gca();   //get axis object
g6.data_bounds(:,1)=[0;fc*2];    //x axis scale
g6.data_bounds(:,2)=[-30;50];    //y axis scale
g6.font_size=5

subplot(4,2,6)
plot(f,20*log10((abs(y3(1:fn)))),"thickness",2)
xgrid()
title("QPSK Spectrum","fontsize",5)
xlabel("Freq [Hz]","fontsize",5)
ylabel("Amplitude [dB]","fontsize",5)
g7=gca();   //get axis object
g7.data_bounds(:,1)=[0;fc*2];    //x axis scale
g7.data_bounds(:,2)=[-30;50];    //y axis scale
g7.font_size=5

subplot(4,2,8)
plot(f,20*log10((abs(y4(1:fn)))),"thickness",2)
xgrid()
title("16QAM Spectrum","fontsize",5)
xlabel("Freq [Hz]","fontsize",5)
ylabel("Amplitude [dB]","fontsize",5)
g8=gca();   //get axis object
g8.data_bounds(:,1)=[0;fc*2];    //x axis scale
g8.data_bounds(:,2)=[-30;50];    //y axis scale
g8.font_size=5
