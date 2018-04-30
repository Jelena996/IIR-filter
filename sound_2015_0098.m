 
% Tacka 1
[x,Fs] = audioread('..\dz2_signali\sound_corrupted.wav');

%crtanje spektrograma
window = 128;
nooverlap = (1/2)*window;
figure;
spectrogram(x, window, nooverlap, [], Fs,'yaxis');
ax.Yscale = 'log';
title('Spektogram originalnog signala');
xlabel('Vreme [s]'),ylabel('Frekvencija [kHz]');
colormap spring;


%Na osnovu prethodno nacrtanog spektrograma vidimo da je potrebno napraviti
%dva noc filtra i jedan NF
%odredjivanje gabarita zeljenih filtara pomocu gotovih funkcija u matlabu

%I) Odredjivanje gabarita za NF filtar
[n,Wp] = ellipord(6000/24000,7000/24000, 0.5,80);
[b,a] = ellip(n,0.5,80,Wp);
y1=filter(b,a,x);
figure;
spectrogram(y1, window, nooverlap, [], Fs,'yaxis');
ax.Yscale = 'log';
title('Spektogram signala posle NF-a');
xlabel('Vreme [s]'),ylabel('Frekvencija [kHz]');
colormap spring;

%II) Odredjivanje gabarita za NO filtre

Wp=[1100 2100 ]/(Fs/2);
Ws=[1400 1900]/(Fs/2);
Rp=0.5;
Rs=50; 
[n,Wp]=ellipord(Wp,Ws, Rp,Rs);
[b,a]=ellip(n,Rp, Rs, Wp,'stop');
y2=filter(b,a,y1);
figure;
spectrogram(y2, window, nooverlap, [], Fs,'yaxis');
ax.Yscale = 'log';
title('Spektogram signala posle prvog noca');
xlabel('Vreme [s]'),ylabel('Frekvencija [kHz]');
colormap spring;

%--------------------------------------------------------------------------

Wp=[2050 3550 ]/(Fs/2);
Ws=[2300 3380]/(Fs/2);
Rp=0.5;
Rs=10;
[n,Wp]=ellipord(Wp,Ws, Rp,Rs);
[b,a]=ellip(n,Rp, Rs, Wp,'stop');

figure;
y3=filter(b,a,y2);
spectrogram(y3, window, nooverlap, [], Fs,'yaxis');
ax.Yscale = 'log';
title('Spektogram signala posle 2. noca');
xlabel('Vreme [s]'),ylabel('Frekvencija [kHz]');
colormap spring;

 
%% Tacke 3,4 i 5

[x,Fs] = audioread('..\dz2_signali\sound_corrupted.wav');

%crtanje spektrograma
window = 128;
nooverlap = (1/2)*window;
figure;
spectrogram(x, window, nooverlap, [], Fs,'yaxis');
ax.Yscale = 'log';
title('Spektogram originalnog signala');
xlabel('Vreme [s]'),ylabel('Frekvencija [kHz]');
colormap spring;


%filtriranje signala pomocu realizovanih funkcija i prethodno odredjenih
%gabarita

[b,a]=design_NF_iit(48000, 6000, 7000, 0.5, 80);
y=filter(b,a,x);

[b,a]=design_BS_bilinear(48000, 2050, 3550, 2300, 3380, 0.5, 70);
y=filter(b,a,y);

[b,a]=design_BS_bilinear(48000, 1100, 2100, 1400, 1900, 0.5, 70);
y=filter(b,a,y);

figure;
spectrogram(y, window, nooverlap, [], Fs,'yaxis');
ax.Yscale = 'log';
title('Spektogram signala posle svih filtara');
xlabel('Vreme [s]'),ylabel('Frekvencija [kHz]');
colormap spring;

filtrirani=audioplayer(10*y,Fs);
filtrirani.play


audiowrite('out_signal_2015_0098.wav', 10*y, Fs);


