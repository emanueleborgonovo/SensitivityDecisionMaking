function [A2,A1,y]=pennz4(x)
%x=[pA5B,pLA,C3,C4];

C1=2;
C2=5;
C3=x(:,3);
C4=x(:,4);
C5=0;
C6=3;
pA5B=x(:,1);
pR5B=0.5/0.17*pA5B;
pTCO=1-pA5B-pR5B;
pLA=x(:,2);
pMA=0.5/0.2*pLA;
pNA=1-pLA-pMA;

if or(pNA<0,pNA>1)
    disp('errorPNA')
end
if or(pTCO<0,pTCO>1)
    disp('errorPTCO')
end

A1=C1;
A2=C2.*pA5B+(C3.*pLA+C4.*pMA+C5.*pNA).*pR5B+pTCO.*max(C3.*pLA+C4.*pMA+C5.*pNA,C6);
y=max(C1,C2.*pA5B+(C3.*pLA+C4.*pMA+C5.*pNA).*pR5B+pTCO.*max(C3.*pLA+C4.*pMA+C5.*pNA,C6));
%y=max(x(:,1),x(:,2).*x(:,7)+(x(:,3).*x(:,10)+x(:,4).*x(:,11)+x(:,5).*x(:,12))*x(:,8)+x(:,9)*max(x(:,3).*x(:,10)+x(:,4).*x(:,11)+x(:,5).*x(:,12),x(:,6)));
end