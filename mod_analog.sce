clear
clc
clf

Fs=10000;   //Sampling frequency
L=10000;    //Length of signal
fb=10; //baseband freq
fc=100;    //carrier freq
m=0.5;      //AM mod. depth
dev=1;     //FM dev. freq

//gen. waveforms
t=0:1/Fs:(L-1)/Fs;  //time vector. 0-1
wave1=sin(2*%pi*fb*t);  //signal
wave2=sin(2*%pi*fc*t);  //carrier
for i=1:L
    am(1,i)=(1+m*wave1(i))*wave2(i);    //AM waveform
end
for i=1:L
    fm(1,i)=sin(2*%pi*fc*t(i)+dev/fb*sin(2*%pi*fb*t(i)));   //FM waveform
end

//fft
f=Fs*(0:(L/2))/L; //freq vector. freq res.=Fs/N*n. n=0-N/2
fn=length(f) //length of freq
win=window('hn',L)      //window func.
y1(1,:)=fft(am.*win);
y2(1,:)=fft(fm.*win);

subplot(2,2,1)
plot(t,am,"thickness",2)
xgrid()
title('AM Waveform',"fontsize",5)
xlabel('Time [s]',"fontsize",5)
ylabel('Amplitude',"fontsize",5)
ga=gca();   //get axis object
ga.data_bounds(:,1)=[0;1/fb*3];    //x axis scale
ga.font_size=5

subplot(2,2,3)
plot(t,fm,"thickness",2)
xgrid()
title('FM Waveform',"fontsize",5)
xlabel('Time [s]',"fontsize",5)
ylabel('Amplitude',"fontsize",5)
gb=gca();   //get axis object
gb.data_bounds(:,1)=[0;1/fb*3];    //x axis scale
gb.font_size=5

//show spectrum
subplot(2,2,2)
plot(f,20*log10(abs(y1(1:fn))),"thickness",2)
xgrid()
title("AM Spectrum","fontsize",5)
xlabel("Freq [Hz]","fontsize",5)
ylabel("Amplitude [dB]","fontsize",5)
gc=gca();   //get axis object
gc.data_bounds(:,1)=[0;fc*2];    //x axis scale
gc.data_bounds(:,2)=[0;100];    //y axis scale
gc.font_size=5

subplot(2,2,4)
plot(f,20*log10((abs(y2(1:fn)))),"thickness",2)
xgrid()
title("FM Spectrum","fontsize",5)
xlabel("Freq [Hz]","fontsize",5)
ylabel("Amplitude [dB]","fontsize",5)
gd=gca();   //get axis object
gd.data_bounds(:,1)=[0;fc*2];    //x axis scale
gd.data_bounds(:,2)=[0;100];    //y axis scale
gd.font_size=5
