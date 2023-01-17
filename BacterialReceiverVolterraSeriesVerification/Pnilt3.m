function fx=Pnilt3(Fx,c);
global N M qd t NT omega
fx(:,:,:,1)=fft(Fx(:,:,:,1),N,3);
fx(:,:,:,2)=N*ifft(Fx(:,:,:,2),N,3);
fx=fx(:,:,1:M,:);
dFx1=size(Fx,1);
dFx2=size(Fx,2);
d=zeros(dFx1,dFx2,qd,2);
e=d;
q=Fx(:,:,N+2:N+qd,:)./Fx(:,:,N+1:N+qd-1,:);
d(:,:,1,:)=Fx(:,:,N+1,:);
d(:,:,2,:)=-q(:,:,1,:);
for r=2:2:qd-1
    w=qd-r;
    e(:,:,1:w,:)=q(:,:,2:w+1,:)-q(:,:,1:w,:)+e(:,:,2:w+1,:);
    d(:,:,r+1,:)=-e(:,:,1,:);
    if r>2
        q(:,:,1:w-1,:)=q(:,:,2:w,:).*e(:,:,2:w,:)./e(:,:,1:w-1,:);
        d(:,:,r,:)=-q(:,:,1,:);
    end
end
A2=zeros(dFx1,dFx2,M,2);
B2=ones(dFx1,dFx2,M,2);
A1=repmat(d(:,:,1,:),[1,1,M,1]);
B1=B2;
z1(1,1,:,1)=exp(-i*omega*t);
z1(1,1,:,2)=conj(z1(1,1,:,1));
z=repmat(z1,[dFx1,dFx2,1]);
for n=2:qd
    Dn=repmat(d(:,:,n,:),[1,1,M,1]);
    A=A1+Dn.*z.*A2;
    B=B1+Dn.*z.*B2;
    A2=A1;
    B2=B1;
    A1=A;
    B1=B;
end
fx=fx+A./B;
fx=sum(fx,4)-repmat(Fx(:,:,1,1),[1,1,M]);
fx=repmat(reshape(exp(c*t)/NT,[1,1,M]),[dFx1,dFx2,1]).*fx;
fx(:,:,1)=2*fx(:,:,1);