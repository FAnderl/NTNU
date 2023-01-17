function ft=NILT1_Volterra(F,tm);
% tm is the time instant (duration);
% F is unknown?  
M=128;
P=3;
Er=1e-8;
alfa=0;
% adjustable parameters
N=2*M;
qd=2*P+1;
NT=2*tm*N/(N-2);
Ken=log(1-1/(1+Er))/NT;
omega=2*pi/NT;
c=alfa-Ken;
Nw=N+qd-1;
Ac=c-1i*omega*(0:Nw);
t=linspace(0,tm,M);
Fq(1,:)=feval(F,Ac);
Fq(2,:)=feval(F,conj(Ac));
ft(1,:)=fft(Fq(1,1:N));
ft(2,:)=N*ifft(Fq(2,1:N));
ft=ft(:,1:M);
d=zeros(2,qd);
e=d;
d(:,1)=Fq(:,N+1);
q=Fq(:,N+2:Nw+1)./Fq(:,N+1:Nw);
d(:,2)=-q(:,1);
for r=2:2:qd-1
    w=qd-r;
    e(:,1:w)=q(:,2:w+1)-q(:,1:w)+e(:,2:w+1);
    d(:,r+1)=-e(:,1);
    if r>2
        q(:,1:w-1)=q(:,2:w).*e(:,2:w)./e(:,1:w-1);
        d(:,r)=-q(:,1);
    end
end
A2=zeros(2,M);
B2=ones(2,M);
A1=repmat(d(:,1),[1,M]);
B1=B2;
z1=exp(-1i*omega*t);
z=[z1;conj(z1)];
for n=2:qd
    Dn=repmat(d(:,n),[1,M]);
    A=A1+Dn.*z.*A2;
    B=B1+Dn.*z.*B2;
    A2=A1;
    B2=B1;
    A1=A;
    B1=B;
end
ft=ft+A./B;
ft=sum(ft)-Fq(2,1);
ft=exp(c*t)/NT.*real(ft);
ft(1)=2*ft(1);
