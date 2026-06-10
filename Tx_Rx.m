fs = 500e3;
PL = 0.005; %[s]
Amp = 0.5;%[V]
t = 0:1/fs:PL - 1/fs;
f = 9e3;
TransSN = 'D25';
RecSC = 'BK81';
savepath = ['C:\Users\pool\Desktop\exp' TransSN '_' RecSC '\']

pulse = Amp*sin(2*pi*f*t);
PRI = 1;
sig = [pulse, zeros(1, (PRI-PL)*fs)]';
repNum = 5;
transSig = repmat(sig, repNum, 1);
filename = [savepath 'f_' num2str(f) 'Hz.mat']
gains.Tx = 20;
gains.TXname = 'B&K 27';
gains.Rx = 30;
gains.RXname = 'B&K 26';

%%
d = daq("ni");
devList = daqlist;
devID = devList.DeviceID(1);

addinput(d, devID, "ai0", "Voltage");
addinput(d, devID, "ai1", "Voltage");

addoutput(d, devID, "ao0", "Voltage");
d.Rate = fs;
preload(d, transSig);
inData = readwrite(d, transSig);

save(filename, "gains", "inData", "transSig");




