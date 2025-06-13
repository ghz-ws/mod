clc
clear
clf

Fs=1000;    //sample freq.
L=1000;     //Length of signal
Baud=100;   //data rate
PN=7        //PN. n=2-20 except n=8,12,13,14,16,19

t=0:1/Fs:(L-1)/Fs;  //time vector. 0-1
f=Fs*(0:(L/2))/L;   //freq vector. freq res.=Fs/L*n. n=0-N/2
fn=length(f) //length of freq
rate=Fs/Baud; //oversample rate

for i=1:L
    if sin(2*%pi*Baud*1/2*t(i))>=0  //rect
        rect(1,i)=1;
    else
        rect(1,i)=0;
    end
end

prbs_lfsr=ones(1,PN);
prbs_tap=[2,3,4,4,6,7,1,6,8,10,1,1,1,15,1,15,12,1,18];  //tap position
for i=1:L/rate
    out=bitxor(prbs_lfsr(PN),prbs_lfsr(prbs_tap(PN-1)-1));
    prbs_lfsr=cat(2,out,prbs_lfsr(1:PN-1));
    for j=1:rate
        if out==0
            prbs(1,rate*(i-1)+j)=0;
        else
            prbs(1,rate*(i-1)+j)=1;
        end
    end
end
if size(prbs)(1)!=L
    prbs(L)=0;
end

//fft
win=window('hn',L)      //window func.
y1(1,:)=fft(rect.*win);
y2(1,:)=fft(prbs.*win);

subplot(2,2,1)
plot(t,rect,"thickness",2)
xgrid()
title("Rect. Waveform","fontsize",5)
xlabel("Time [s]","fontsize",5)
ylabel("Amplitude","fontsize",5)
g1=gca();   //get axis object
g1.data_bounds(:,1)=[0;20/Baud];   //y axis scale
g1.data_bounds(:,2)=[0;1.2];   //y axis scale
g1.font_size=5

subplot(2,2,3)
plot(t,prbs,"thickness",2)
xgrid()
title("PRBS Waveform","fontsize",5)
xlabel("Time [s]","fontsize",5)
ylabel("Amplitude","fontsize",5)
g2=gca();   //get axis object
g2.data_bounds(:,1)=[0;20/Baud];   //y axis scale
g2.data_bounds(:,2)=[0;1.2];   //y axis scale
g2.font_size=5

subplot(2,2,2)
plot(f,20*log10(abs(y1(1:fn))),"thickness",2)
xgrid()
title("Rect. Spectrum","fontsize",5)
xlabel("Freq [Hz]","fontsize",5)
ylabel("Amplitude [dB]","fontsize",5)
g3=gca();   //get axis object
g3.data_bounds(:,2)=[0;60];    //y axis scale
g3.font_size=5

subplot(2,2,4)
plot(f,20*log10(abs(y2(1:fn))),"thickness",2)
xgrid()
title("PRBS Spectrum","fontsize",5)
xlabel("Freq [Hz]","fontsize",5)
ylabel("Amplitude [dB]","fontsize",5)
g4=gca();   //get axis object
g4.data_bounds(:,2)=[0;60];    //y axis scale
g4.font_size=5
