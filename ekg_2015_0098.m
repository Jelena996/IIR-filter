
load('..\dz2_signali\ecg_corrupted');
x=val;
Fs=360;
dt=1/Fs;
Tuk=length(x)/Fs;
t=0:dt:Tuk-dt;
figure
plot(t,x), title('Originalni signal prije filtriranja');
xlabel('t'), ylabel('x');

%Zadatak 3
%filtriranje signala projektovanim VF filtrom

[b,a]=baseline_drift_filter(360,0.4,1,30,0.5);
x=filter(b,a,x);
figure
plot(t,x), title('Signal nakon filtriranja VF filtrom');
xlabel('t'), ylabel('x');


%Zadatak 4
%filtriranje signala projektovanim NO filtrom
[b,a]=power_line_noise_filter(360, 60, 40, 0.5);
x=filter(b,a,x);
figure
plot(t,x), title('Signal nakon filtriranja noc filtrom');
xlabel('t'), ylabel('x');

%Zadatak 5
%Amplituske i fazne k-ke se crtaju pozivanjem prethodih funkcija






