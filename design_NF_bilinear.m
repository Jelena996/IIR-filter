function [ bd, ad ] = design_NF_bilinear( Fs, Fp, Fa, Ap, Aa )
%projektuje NF filtar koriscenjem eliptickog analognog prototipa i
%bilinearne transformacije
Nfreqz=10000;
Wa=2*pi*Fa/Fs;
Wp=2*pi*Fp/Fs;


%gabariti normalizovanog analognog filtra
%1) predistorzija

WaPred=2*Fs*tan(Wa/2);
WpPred=2*Fs*tan(Wp/2);

wpN=1;
waN=WaPred/WpPred; 
Aav=Aa;
Apv=Ap;

while(1)
   
    k=sqrt(1-(wpN/waN)^2);
    D=(10^(0.1*Aav)-1)/(10^(0.1*Apv)-1);
    q0=(1/2)*((1-sqrt(k))/(1+sqrt(k)));
    q=q0+2*q0^5+15*q0^9+15*q0^13;
    N=ceil(log10(16*D)/(log10(1/q)));    
    
    [z,p,k]=ellipap(N,Apv,Aav);
    baN=k*poly(z);aaN=poly(p);
    
%transformacija NF -> NF
    [b,a] = lp2lp(baN,aaN,WpPred);
    
%diskretizacija
    [bd,ad]=bilinear(b,a,Fs);
    df=(Fs/2)/Nfreqz;
    [hd,wd]= freqz(bd,ad,Nfreqz);
    Hd=abs(hd);
%prvojera gabarita
    NPOK = 0;
    POOK = 0;
     
   ia=floor(Fa/df)+1;
   ip=ceil(Fp/df)+1;
   Hdprov=Hd(1:ip);
   Haprov=Hd(ia:length(Hd));

    if(max(20*log10(Haprov))>-Aa )
        Aav=Aav+0.1;
    else
        NPOK=1;
    end
    
    if(min(20*log10(Hdprov))<-Ap)
        Apv=Apv-0.1;
    else
        POOK=1;
    end
    
    if((NPOK==1)&&(POOK==1))
         N
        break
    end

end
  
%crtanje amplitudske karakteristike NF filtra
figure;
f=(wd*Fs)/(2*pi);
plot(f,20*log10(Hd)),xlabel('Ucestanost (Hz)'), ylabel('20log|H|','LineWidth', 2),title('Amplitudksa k-ka NF filtra');


%crtanje gabarita NF filtra
hold on
xh = [Fp/10 Fp]; yh = [-Ap -Ap];
xv = [Fp Fp]; yv = [-Ap 0];
x2h = [Fa Fa*10]; y2h = [-Aa -Aa];
x2v = [Fa Fa]; y2v = [-Aa -2*Aa];
plot(xh,yh,'r',xv,yv,'r',x2h,y2h,'r',x2v,y2v,'r');
hold off


%crtanje fazne k-ke NF filtra
figure
faza=phase(hd');
plot(f,faza), title('Fazna k-ka NF filtra');
xlabel('Ucestanost (Hz)');ylabel('phase [rad]');















end

