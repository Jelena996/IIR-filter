function [ bd,ad ] =design_BS_bilinear( fs,fp1,fp2,fa1,fa2,Ap,Aa )
%projektuje BS filtar koriscenjem eliptickog analognog prototipa i
%bilinearne transformacije

wp1=2*pi*fp1;
wp2=2*pi*fp2;
wa1=2*pi*fa1;
wa2=2*pi*fa2;
Aav=Aa;
Apv=Ap;
T=1/fs;


%odredjivanje specifikacija digitalnog PO u slucaju da uneseni parametri ne
%zadovoljavaju uslov geometrijske simetrije granicnih ucestanosti
w0=sqrt(wp1*wp2);
if (wa2 > w0^2/wa1)
   wa2 = w0^2/wa1; 
   fa2=wa2/(2*pi);
end

if(wa1 < w0^2/wa2)
 wa1 = w0^2/wa2; 
 fa1=wa1/(2*pi);
end


%odredjivanje specifikacija normalizovanog analognog prototipa
wa1Pred = 2/T*tan(fa1*2*pi/fs/2);
wa2Pred = 2/T*tan(fa2*2*pi/fs/2);
wp1Pred = 2/T*tan(fp1*2*pi/fs/2);
wp2Pred = 2/T*tan(fp2*2*pi/fs/2);

%kako sad imamo nove parametre za analogni filtar zbog predistorzije,
%moramo izracunati i novo B i w0
B = wp2Pred-wp1Pred;
w0 = sqrt(wp1Pred*wp2Pred);

%uzimanje strozijeg uslova kako bi bio zadovoljen usloc
%w0^2=sqrt(wa1Prd*wa2Pred); tako da jednu ucestanost zadrzavamo, a drugu
%proracunavamo

if (wa1Pred > w0^2/wa2Pred)
   wa1Pred = w0^2/wa2Pred; 
end

%gabariti normalizovanog prototipa
%NF prototip
wpN=1;
waN=wa1Pred*B/(w0^2-wa1Pred^2);

I=0;
while(1)
    %parametri za odredjivanje analognog prototipa preko
    %elipticke aproksimacije(ova aproksimacija je izabrana jer daje
    %najmanji red filtra, iako u zadatku nije trazeno
k=sqrt(1-(wpN/waN)^2);
    D=(10^(0.1*Aav)-1)/(10^(0.1*Apv)-1);
    q0=(1/2)*((1-sqrt(k))/(1+sqrt(k)));
    q=q0+2*q0^5+15*q0^9+15*q0^13;
    N=ceil(log10(16*D)/(log10(1/q)));    
    I=I+1;
    I
N
[z,p,k]=ellipap(N,Apv,Aav);
baN=k*poly(z);aaN=poly(p);

%Transformacija normalizovanog prototipa u NO
[ba,aa]=lp2bs(baN,aaN,w0,B);
n=numel(aa);

%diskretizacija    
[bd,ad]=bilinear(ba,aa,fs);
Nfreqz = 10000;
[hd,Wd]=freqz(bd,ad,Nfreqz);
Hd=abs(hd);

%provjera gabarita
 NPOK = 0;
 POOK = 0;
 df=(fs/2)/Nfreqz;
 ia1=floor(fa1/df)+1;
 ia2=ceil(fa2/df)+1;
 ip1=ceil(fp1/df)+1;
 ip2=floor(fp2/df)+1;

 Ha=Hd(ia1:ia2);   %Amplitudska Karakteristika Nepropusnog Opsega
 Hp=[Hd(1:ip1)' Hd(ip2:length(Hd))'];   %Amplitudska Karakteristika Propusnog Opsega

    if(max(20*log10(Ha))>-Aa)         
        Aav=Aav+0.1;
    else
        NPOK=1;
    end
    
    if(min(20*log10(Hp))<-Ap)
        Apv=Apv-0.1;
    else
        POOK=1;
    end   
    
    if((NPOK==1)&&(POOK==1))
        n=filtord(bd,ad);
        n %dobro je, paran je
        N
        break
    end
 
end


%amplitudska karakteristika digitalnog filtra u logaritamskoj razmeri
f = Wd/(2*pi)*fs;
figure;
plot(f,20*log10(Hd), 'LineWidth', 2);
title('Amplitudska karakteristika BS filtra');
xlabel('Ucestanost (Hz)');
ylabel('20log|H|');

%crtanje gabarita
hold on
xh = [fa1 fa2]; yh = [-Aa -Aa];
xv = [fa1 fa1]; yv = [-Aa -2*Aa];
x2v = [fa2 fa2]; y2v = [-Aa -2*Aa];

x2h = [fp1 fp1/10]; y2h = [-Ap -Ap];
x3v = [fp1 fp1]; y3v = [-Ap 0];

x3h = [fp2 fp2*10]; y3h = [-Ap -Ap];
x4v = [fp2 fp2]; y4v = [-Ap 0];

plot(xh,yh,'r',xv,yv,'r',x2h,y2h,'r',x2v,y2v,'r', x3v,y3v,'r',x3h,y3h,'r',x4v,y4v,'r')
hold off

%crtanje fazne k-ke
figure
faza=phase(hd');
plot(f,faza), title('Fazna k-ka BS filtra');
xlabel('Ucestanost (Hz)');ylabel('phase [rad]');

end

