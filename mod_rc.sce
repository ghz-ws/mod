clear
clc
clf

Fs=1000     //sample freq
Baud=100    //symbol rate
alpha=0.5   //roll off
TAP=81      //TAP
L=2000;     //Length of signal
PN=9        //PRBS n. n=2-20 except n=8,12,13,14,16,19
D=2         //bits/symbol. 1,2,4,6

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

[hm,fr]=frmag(coe,200); //calc freq response
fr_nom=Fs*fr;

subplot(2,3,4)
plot(fr_nom,hm,"thickness",2)  //plot freq response
xgrid()
title("RC Freq Response","fontsize",5)
xlabel("Symbol Rate [sps]","fontsize",5)
ylabel("Amplitude","fontsize",5)
cp=gca();   //get axis object
cp.data_bounds(:,1)=[0;Baud];    ///x axis scale
cp.font_size=5

//gen. time vector and baseband
for i=1:L
    X(i)=(i-1)*1/Fs;    //time vector
end
prbs_lfsr=ones(1,PN);
prbs_tap=[2,3,4,4,6,7,1,6,8,10,1,1,1,15,1,15,12,1,18];  //tap position
rate=Fs/Baud;   //samples per symbol
for i=1:L/rate
    for k=1:D
        seed=bitxor(prbs_lfsr(PN),prbs_lfsr(prbs_tap(PN-1)-1));
        prbs_lfsr=cat(2,seed,prbs_lfsr(1:PN-1));
        temp(k)=seed;
    end
    for j=1:rate
        if D==1 //1bit/symbol ASK/PSK
            x1(rate*(i-1)+j)=(temp(1)-0.5)*2;
            if j==1
                I(rate*(i-1)+j)=(temp(1)-0.5)*2;   //format 1/-1
                Q(rate*(i-1)+j)=0;
            else
                I(rate*(i-1)+j)=0;
                Q(rate*(i-1)+j)=0;
            end
        elseif D==2 //2bit/symbol QPSK
            x1(rate*(i-1)+j)=(temp(1)-0.5)*2;
            if j==1
                I(rate*(i-1)+j)=(temp(1)-0.5)*2;   //format 1/-1
                Q(rate*(i-1)+j)=(temp(2)-0.5)*2;
            else
                I(rate*(i-1)+j)=0;
                Q(rate*(i-1)+j)=0;
            end
        elseif D==4 //4bit/symbol 16QAM
            i_base=temp(1)+2*temp(2);   //0-3
            q_base=temp(3)+2*temp(4);
            x1(rate*(i-1)+j)=(i_base-1.5)/1.5;
            if j==1
                I(rate*(i-1)+j)=(i_base-1.5)/1.5;   //format 1/0.3/-0.3/-1
                Q(rate*(i-1)+j)=(q_base-1.5)/1.5;
            else
                I(rate*(i-1)+j)=0;
                Q(rate*(i-1)+j)=0;
            end
        elseif D==6 //6bit/symbol 64QAM
            i_base=temp(1)+2*temp(2)+4*temp(3);   //0-7
            q_base=temp(4)+2*temp(5)+4*temp(6);
            x1(rate*(i-1)+j)=(i_base-1.5)/1.5;
            if j==1
                I(rate*(i-1)+j)=(i_base-3.5)/3.5;
                Q(rate*(i-1)+j)=(q_base-3.5)/3.5;
            else
                I(rate*(i-1)+j)=0;
                Q(rate*(i-1)+j)=0;
            end
        end 
    end
end

//process filter
lfsr=zeros(1,size(coe)(1));
for i=1:L
    lfsr=cat(2,I(i),lfsr(1:size(coe)(1)-1));
    I_OUT(i)=sum(coe'.*lfsr);
end
lfsr=zeros(1,size(coe)(1));
for i=1:L
    lfsr=cat(2,Q(i),lfsr(1:size(coe)(1)-1));
    Q_OUT(i)=sum(coe'.*lfsr);
end

subplot(2,2,2)
plot(X(1:L-TAP/2),I(1:L-TAP/2),X(1:L-TAP/2),I_OUT(TAP/2+1:L),"thickness",2)
legend("Input", "Output", 1)
xgrid()
title("Waveform","fontsize",5)
xlabel("Time [s]","fontsize",5)
ylabel("Amplitude","fontsize",5)
gt=gca();   //get axis object
gt.data_bounds(:,1)=[0;10/Baud];    //x axis scale
gt.font_size=5
e2=gce();
e2.font_size=5

//fft
f=Fs*(0:(L/2))/L; //freq vector. freq res.=Fs/N*n. n=0-N/2
fn=length(f) //length of freq
win=window('hn',L)      //window func.
y1(1,:)=fft(x1'.*win);
y2(1,:)=fft(I_OUT'.*win);

subplot(2,2,1)
plot(f,20*log10(abs(y1(1:fn))),f,20*log10(abs(y2(1:fn))),"thickness",2)
xgrid()
legend("Input", "Output", 1)
title("Spectrum","fontsize",5)
xlabel("Freq [Hz]","fontsize",5)
ylabel("Amplitude [dB]","fontsize",5)
gs=gca();   //get axis object
gs.data_bounds(:,1)=[0;Baud*4];    //x axis scale
gs.data_bounds(:,2)=[-40;50];    //y axis scale
gs.font_size=5
es=gce();
es.font_size=5

subplot(2,3,5)
tv=1:1:rate*2+1;    //time vector for eye pattern
for i=1:(L-TAP)/rate/2
    plot(tv',I_OUT(rate*2*(i-1)+TAP:rate*2*(i-1)+TAP+rate*2),"thickness",2)
end
xgrid()
title("Eye pattern","fontsize",5)
xlabel("Time [sample]","fontsize",5)
ylabel("Amplitude","fontsize",5)
et=gca();   //get axis object
et.data_bounds(:,1)=[1;rate*2+1];    //x axis scale
et.font_size=5

subplot(2,3,6)
plot(I_OUT(TAP/2+1:L),Q_OUT(TAP/2+1:L),"thickness",2)
xgrid()
title("Constellation","fontsize",5)
xlabel("I","fontsize",5)
ylabel("Q","fontsize",5)
eiq=gca();   //get axis object
eiq.font_size=5
