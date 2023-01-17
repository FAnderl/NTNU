function [ft,t]=NILT3_Volterra(F,tm);
global N M qd t NT omega
M=128;
P=3;
Er=1e-8;
dim=3;
alfa=zeros(dim,1); % adjustable parameters
N=2*M;
qd=2*P+1;
NT=2*tm*N/(N-2);
Ken=log(1-1/nthroot(1+Er,3))/NT;
omega=2*pi/NT;
c=alfa-Ken;
Nw=N+qd;
Ac=repmat(c,1,Nw)-repmat(1i*omega*(0:Nw-1),dim,1);
Acc=conj(Ac);
roz=1:Nw;
rozc=Nw+1:2*Nw;
t=linspace(0,tm,M);
Fs=zeros(2*Nw,2*Nw,Nw,2);
[s1,s2,s3]=ndgrid([Ac(1,:),Acc(1,:)],[Ac(2,:),Acc(2,:)],Ac(3,:));
Fs(:,:,:,1)=feval(F,s1,s2,s3);
[s1,s2,s3]=ndgrid([Ac(1,:),Acc(1,:)],[Ac(2,:),Acc(2,:)],Acc(3,:));
Fs(:,:,:,2)=feval(F,s1,s2,s3);
Fs=Pnilt3(Fs,c(3));
Fst=zeros(2*Nw,Nw,M,2);
Fst(:,:,:,1)=Fs(:,roz,:);
Fst(:,:,:,2)=Fs(:,rozc,:);
Fst=permute(Fst,[1,3,2,4]);
Fs=Pnilt3(Fst,c(2));
Fs=ipermute(Fs,[1,3,2]);
Fst=zeros(Nw,M,M,2);
Fst(:,:,:,1)=Fs(roz,:,:);
Fst(:,:,:,2)=Fs(rozc,:,:);
Fst=permute(Fst,[3,2,1,4]);
Fs=Pnilt3(Fst,c(1));
Fs=ipermute(Fs,[3,2,1]);
ft=zeros(1,M);
for k=1:M
    ft(k)=real(Fs(k,k,k));
end