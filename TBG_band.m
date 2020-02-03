clear;monitor=0;
%G vectors
cutoff=4;waven=(2*cutoff+1)^2;
G=zeros(waven,2);
m=1;
for i=1:2*cutoff+1
  for j=1:2*cutoff+1
     G(m,:)=[i-cutoff-1 j-cutoff-1];
     m=m+1;
  end
end
%
valley=1;                      %valley index
a0=2.46;                       %lattice constant
vf=2135.4*a0;u1=79.7;u2=97.5;  %fermi velocity and interlayer tunnelling
Delta=0;                       %stagger potential by hBN
theta=1.05*pi/180;             %twist angle                                                                                                                                                   
sx=[0,1;1,0];sy=[0,-1i;1i,0];sz=[1,0;0,-1];s0=[1,0;0,1]; %pauli matrix
a1=[1;0];a2=[1/2;sqrt(3)/2];                            
b=4*pi/(3*a0)*[1;0];
Rot1=[cos(theta/2),sin(theta/2);-sin(theta/2),cos(theta/2)];
Rot2=[cos(theta/2),-sin(theta/2);sin(theta/2),cos(theta/2)];
%tunnelling matrix
fai=valley*2*pi/3;
T1=[u1,u2;u2,u1];
T2=[u1,u2*exp(-1i*fai);u2*exp(1i*fai),u1];
T3=[u1,u2*exp(1i*fai);u2*exp(-1i*fai),u1];
%moire vector
bm=norm(Rot1*b-Rot2*b);
G1M=[-1/2;-sqrt(3)/2]*sqrt(3)*bm;
G2M=[1;0]*sqrt(3)*bm;
K1=valley*[sqrt(3)/2;-1/2]*bm;
K2=valley*[sqrt(3)/2;1/2]*bm;
%k path
i=1;
for m=0:1:100  k(:,i)=[-sqrt(3)*(m-100)/200 -(m-100)/200] ;i=i+1; end
for m=101:1:200 k(:,i)=[(m-100)*sqrt(3)/200 0] ; i=i+1;end
for m=201:1:300 k(:,i)=[sqrt(3)/2 -(m-200)/200] ; i=i+1;end
%
Energy=zeros(waven*4,301);
for d=1:301
 Hint=zeros(waven*2,waven*2);
 H1=zeros(waven*2,waven*2);
 H2=zeros(waven*2,waven*2);
 km=k(:,d)*bm;
 for i=1:waven
   g=G(i,:);
   kvec1=Rot2*(km-K1+g(1)*G1M+g(2)*G2M);kvec2=Rot1*(km-K2+g(1)*G1M+g(2)*G2M);
   H1(2*i-1:2*i,2*i-1:2*i)=-vf*(kvec1(1)*valley*sx+kvec1(2)*sy)+Delta*sz;
   H2(2*i-1:2*i,2*i-1:2*i)=-vf*(kvec2(1)*valley*sx+kvec2(2)*sy); 
   for j=1:waven
      g1=G(i,:);g2=G(j,:);
      if i==j
        Hint(2*i-1:2*i,2*i-1:2*i)=T1;
      end
      if g2(1)-g1(1)==-valley*1 && g1(2)==g2(2)
        Hint(2*i-1:2*i,2*j-1:2*j)=T2; 
      end
      if g2(1)-g1(1)==-valley*1 && g2(2)-g1(2)==-valley*1
        Hint(2*i-1:2*i,2*j-1:2*j)=T3; 
      end     
  end
 end

 H=[H1 Hint;Hint' H2];
 e=eig(H);
 Energy(:,d)=e;
 monitor=monitor+1
end

for i=waven*2-10:waven*2+10
    hold on;
    plot(Energy(i,:),'-','Linewidth',1.1,'Color','b');
    ylabel('E(meV)');
    set(gca,'xticklabel',[])
    box on
end
axis([0 301 -80 100]);
set(gca,'XTick',0:100:300);
set(gca,'YTick',-80:20:100);
set(gca,'XTickLabel',{'K_+^m';'\Gamma^m';'M^m';'K_-^m'});
hold on 
plot([100 100],[-200 200],'--','LineWidth',0.9,'Color','k');
hold on
plot([200 200],[-200 200],'--','LineWidth',0.9,'Color','k');
