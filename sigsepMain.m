%{
This is a template main script for signal separation problem 

%}
clc; clear all; close all;

[music,Fs1] = audioread('musicf1.wav');
%sound(musicw,Fs1);
[speech,Fs2] = audioread('speechf1.wav');
%sound(speechw,Fs2);
[mixed,Fs3] = audioread('mixedf1.wav');
%sound(mixed,Fs3);   % Audio files are read seperately.

% ADD CODE TO COMPUTE MAG. SPECTOROGRAMS OF MUSIC AND SPEECH
music_spectrum=stft(music', 2048, 256, 0, hann(2048));
MagMusic=abs(music_spectrum);
PhaseMusic=music_spectrum./MagMusic;

speech_spectrum=stft(speech', 2048, 256, 0, hann(2048));
MagSpeech=abs(speech_spectrum);
PhaseSpeech=speech_spectrum./MagSpeech;     

% ADD CODE TO CPMPUTE MAG. SPECTROGRAM AND PHASE OF MIXED
mixed_spectrum=stft(mixed', 2048, 256, 0, hann(2048));
MagMixed=abs(mixed_spectrum);
PhaseMixed=mixed_spectrum./MagMixed;    % Magnitude spectrograms and Phases are computed seperately.

niter = 250;

Bminit = load('Bminit.mat','Bm'); Bminit=Bminit.Bm;
Wminit = load('Wminit.mat','Wm'); Wminit=Wminit.Wm;
Bsinit = load('Bsinit.mat','Bs'); Bsinit=Bsinit.Bs;
Wsinit = load('Wsinit.mat','Ws'); Wsinit=Wsinit.Ws; 
%Initial base and weight matrices are load seperately.


[Bm,Wm,error1] = doNMF(MagMusic,niter,Bminit,Wminit); 
figure; plot(1:niter,error1);
xlabel('iteration number'); ylabel('Mean Squared Error'); title('Music signal');

[Bs,Ws,error2] = doNMF(MagSpeech,niter,Bsinit,Wsinit);
figure; plot(1:niter,error2);
xlabel('iteration number'); ylabel('Mean Squared Error'); title('Speech signal');
%With doNMF fuction, Optimal bases and weights are computed.

save("BasesForMusic(NMF).mat","Bm"); save("BasesForSpeech(NMF).mat","Bs");

% ADD CODE TO SEPARATE SIGNALS
[music_recv, speech_recv] = separate_signals(MagMixed,Bm,Bs,niter);
% signals are seperated from the mixed audio file with the calculated optimum base matrices(Bm,Bs).
 

% % ADD CODE TO MULTIPLY BY PHASE AND RECONSTRUCT TIME DOMAIN SIGNAL 
seperatedMusic=stft((music_recv.*PhaseMixed),2048,256,0,hann(2048));
seperatedMusic=transpose(seperatedMusic);
sound(seperatedMusic,Fs3);
%
pause(20)   %so that the voices do not interfere.
%  
seperatedSpeech=stft((speech_recv.*PhaseMixed),2048,256,0,hann(2048));
seperatedSpeech=transpose(seperatedSpeech);
sound(seperatedSpeech,Fs3);

%Seperated signals are taken to the time domain(Inverse Stft). Then, 
%Then, they are played so that we can observe if the algorithm is working well.

% WRITE TIME DOMAIN SPEECH AND MUSIC USING audiowrite with 16000 sampling
Fs=16000;
audiowrite('seperatedMusic.wav',seperatedMusic,Fs);
audiowrite('seperatedSpeech.wav',seperatedSpeech,Fs); % Finally, they are saved as waw files with 16k samples.





