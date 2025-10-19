C1=2;
C2=5;
C3=10.3;
C3min=5.15;
C3max=15.45;
C4=5;
C4min=2.5;
C4max=7.5;
C5=0;
C6=3;
pA5B=0.17;
pA5Bmin=0.09;
pA5Bmax=0.255;
pR5B=0.5/017*pA5B;
pTCO=1-pA5B-pR5B;
pLA=0.2;
pLAmin=0.1;
pLAmax=0.3;
pMA=0.5/0.2*pLA;
pNA=1-pLA-pMA;
x=[C1,C2,C3,C4,C5,C6,pA5B,pR5B,pTCO,pLA,pMA,pNA];

xhigh=[C1,C2,C3max,C4max,C5,C6,pA5Bmax,pR5B,pTCO,pLAmax,pMA,pNA];
xlow=[C1,C2,C3min,C4min,C5,C6,pA5Bmin,pR5B,pTCO,pLAmin,pMA,pNA]

% For checking: yy=max(C1,C2*pA5B+(C3*pLA+C4*pMA+C5*pNA)*pR5B+pTCO*max((C3*pLA+C4*pMA+C5*pNA),C6))

[A1,A2,y]=pennz(x)
[k,U,DX,y,ff,phi]=finitechangespennz(x,xhigh);


