function [k,U,DX,y,ff,phi]=finitechangespennz4(x,x1)
%x0= base case; x1=nuovo punto
dx=x1-x;
n=size(x,2);
m=sumbincoeff(n);
U=zeros(n,m);
for i=1:n
    for j=1:m
        if i==j
            U(i,j)=1;
        end
    end 
end
k=0;
k=n;
for l=2:n
    a=combnk(1:n,l);
    b=binomial(n,l);
    counter=k;
    for k=k+1:k+b;
        if k==m-1
            U(:,k)=1;
        elseif k==m
            U(:,k)=0;
        else        
            for z=1:n
                for w=1:size(a,2)
                    if z==a(k-counter,w)
                    U(z,k)=1;
                    end
                end
            end
        end
    end
end
U;
%Base case value is in the last column of U

% Creazione della matrice dU
for i=1:k+1
    for j=1:n
        dU(j,i)=dx(j)*U(j,i);
    end
end
% Creazione della matrice DX
for i=1:k+1
    for j=1:n
        DX(j,i)=dx(j)*U(j,i)+x(j);
    end
end

y=ones(1,m);
%Valutazione del modello esterno
for i=1:m
    y(i)=pennz4(DX(:,i)');
end
y;
% Ortogonalizzazione fij
i=0;
ff=zeros(1,m);
count=0;
k=0;
for l=1:n
    a=combnk(1:n,l);
    b=binomial(n,l);
    counter=k; 
    for i=k+1:k+b;
        ff(i)=y(i);
        for u=1:k
            ep=0;
            for j=1:n
                ep=U(j,i)*U(j,u)+ep;
            end
            %ep=ep
            %i=i
            %u=u
            %pause
            if ep==sum(U(:,u))
                ind=1;
            else
                ind=0;
            end
            ff(i)=ff(i)-ind*ff(u);
        end
    ff(i)=ff(i)-y(m);
    end
    k=k+b;
end
ff;
%Importance
phi=zeros(1,n);
for j=1:n
    for i=1:m
        phi(j)=phi(j)+ff(i)*U(j,i);
    end
end
phi;
%Importance normalized

