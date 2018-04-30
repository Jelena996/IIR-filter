function [bd, ad] =  baseline_drift_filter(fs, fa, fp, Aa, Ap )
%projektuje IIR filtar propusnik visokih uèestanosti korišæenjem normalizovanog Batervortovog prototipa. 
%Funkcija kao argumente prima uèestanost odabiranja fs, graniènu uèestanost nepropusnog opsega fa, 
%graniènu uèestanost nepropusnog opsega fp i odgovarajuæa slabljenja u nepropusnom (?a) i propusnom (?p) opsegu.
%Kao povratnu vrednost, funkcija vraæa koeficijente polinoma u imeniocu I brojiocu funkcije prenosa 
%dobijenog filtra.

Nfreqz=180000;
wp=2*pi*fp/fs;
wa=2*pi*fa/fs;

%parametri analognog prototipa sa predistorzijom
Wp=2*fs*tan(wp/2);
Wa=2*fs*tan(wa/2);
  
%parametri za normalizovani filtar 
Wpp=1;    
Apnovo=Ap;
Aanovo=Aa;
Wap=Wp/Wa;

while(1)
epsilon=sqrt(10^(0.1*Apnovo)-1);
D = (10^(0.1*Aanovo)-1)/(10^(0.1*Apnovo)-1);
k = Wpp/Wap;
N = ceil( log10(D)/2/log10(1/k) );

skalFaktorPolova=nthroot(epsilon,N);
[z,p,k]=buttap(N);
p = p/skalFaktorPolova;
k=k/epsilon;
b_prot=k*poly(z);
a_prot=poly(p);

%transformacija NF -> VF
[b,a] = lp2hp(b_prot,a_prot,Wp);

%diskretizacija
[bd,ad]=bilinear(b,a,fs);
[hd,wd]= freqz(bd,ad,Nfreqz);
Hd=abs(hd);

df=(fs/2)/Nfreqz;
%prvojera gabarita

    NPOK = 0;
    POOK = 0;
     
   ia=ceil(fa/df)+1;
   ip=floor(fp/df)+1;
   Hdprov=Hd(ip:length(Hd));
   Haprov=Hd(1:ia);

    if(max(20*log10(Haprov))>-Aa )
        Aanovo=Aanovo+0.1;
    else
        NPOK=1;
    end
    if(min(20*log10(Hdprov))<-Ap)
        Apnovo=Apnovo-0.1;
    else
        POOK=1;
    end   
    if((NPOK==1)&&(POOK==1))
        N
        break
    end

end
  
%crtanje amplitudske karakteristike VF filtra
figure
f=(wd*fs)/(2*pi);
plot(f,20*log10(Hd)),xlabel('Ucestanost (Hz)'), ylabel('20log|H|','LineWidth', 2), title('Amplitudska k-ka VF filtra u logaritamskoj razmjeri');

%crtanje gabarita
hold on
xh = [fp fs/2]; yh = [-Ap -Ap];
xv = [fp fp]; yv = [-Ap 0];
x2h = [fa fa/10]; y2h = [-Aa -Aa];
x2v = [fa fa]; y2v = [-Aa -2*Aa];
plot(xh,yh,'r',xv,yv,'r',x2h,y2h,'r',x2v,y2v,'r')
hold off

%crtanje fazne k-ke
figure
faza=phase(hd');
plot(f,faza), title('Fazna k-ka VF filtra');
xlabel('Ucestanost (Hz)');ylabel('phase [rad]');

end





