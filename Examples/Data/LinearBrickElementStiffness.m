
%this files is used to verify the result of TestStiffness.cpp
Ke=Stiffnesske(210E9,0.3,0,0,0,1,0,0,1,1,0,0,1,0,0,0,1,1,0,1,1,1,1,0,1,1)
Ke=Stiffnesske(210E9,0.3,0,0,0,1,0,0,1,0.8,0,0,1,0,0,0,1,1,0,1,0.8,0.8,0.8,0,1,1)

function Ke=Stiffnesske(E,NU,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8)
%��Ԫ�նȾ�����㺯�������ݸ����ĸ���������㵥Ԫ�նȾ���
%����E,NUΪ����ģ���Ͳ��ɱȣ�x1-z8Ϊ��Ԫ�˽ڵ�����꣬���KeΪ��Ԫ�նȾ���
%**************��д��������  ������ͨ��ѧ ��еѧԺ ˶1005***************

Loc=[x1 y1 z1;x2 y2 z2;x3 y3 z3;x4 y4 z4;x5 y5 z5;x6 y6 z6;x7 y7 z7;x8 y8 z8;];

gsx=[-0.7745966692 0 0.7745966692]; %��˹����������ϵ��
gsw=[0.55555555556 0.888888888889 0.55555555556];
Ke=zeros(24,24);
for ii=1:3 %��ά��˹���
sx=gsx(ii);
sw=gsw(ii);
for jj=1:3
nx=gsx(jj);
nw=gsw(jj);
for kk=1:3
tx=gsx(kk);
tw=gsw(kk);
Ke=Ke+sw*nw*tw*BDcalc(sx,nx,tx,Loc,E,NU);
end
end
end
end

function BD=BDcalc(s,n,t,Loc,E,NU)
%���ݸ�����s,n,t��ֵ���㺯��ֵBD��������Χ�Ļ���
%**************��д��������  ������ͨ��ѧ ��еѧԺ ˶1005***************

N1=(1-s)*(1-n)*(1-t)/8;
N2=(1+s)*(1-n)*(1-t)/8;
N3=(1+s)*(1+n)*(1-t)/8;
N4=(1-s)*(1+n)*(1-t)/8;
N5=(1-s)*(1-n)*(1+t)/8;
N6=(1+s)*(1-n)*(1+t)/8;
N7=(1+s)*(1+n)*(1+t)/8;
N8=(1-s)*(1+n)*(1+t)/8;

dNsnt=[-(1-n)*(1-t)/8,  -(1-s)*(1-t)/8,  -(1-s)*(1-n)/8; %N1-8��s,n,t�ĵ�������
(1-n)*(1-t)/8,   -(1+s)*(1-t)/8,  -(1+s)*(1-n)/8;
(1+n)*(1-t)/8,   (1+s)*(1-t)/8,  -(1+s)*(1+n)/8;
-(1+n)*(1-t)/8,  (1-s)*(1-t)/8,  -(1-s)*(1+n)/8;
-(1-n)*(1+t)/8,  -(1-s)*(1+t)/8,   (1-s)*(1-n)/8;
(1-n)*(1+t)/8,   -(1+s)*(1+t)/8,   (1+s)*(1-n)/8;
(1+n)*(1+t)/8,   (1+s)*(1+t)/8,   (1+s)*(1+n)/8;
-(1+n)*(1+t)/8,  (1-s)*(1+t)/8,   (1-s)*(1+n)/8;];
dNsnt=dNsnt';
J=dNsnt*Loc;
detJ=det(J);

dNxyz=J\dNsnt;
B=zeros(6,24);
for ii=1:8 %����B����
Bii=[dNxyz(1,ii) 0 0;0 dNxyz(2,ii) 0;0 0 dNxyz(3,ii);
dNxyz(2,ii) dNxyz(1,ii) 0;
0 dNxyz(3,ii) dNxyz(2,ii);
dNxyz(3,ii) 0 dNxyz(1,ii);];
B(:,3*(ii-1)+1:3*ii)=Bii;
end

D=[1-NU NU NU 0 0 0;NU 1-NU NU 0 0 0;NU NU 1-NU 0 0 0;0 0 0 0.5-NU 0 0;0 0 0 0 0.5-NU 0;0 0 0 0 0 0.5-NU;];
D=D*(E/((1+NU)*(1-2*NU))); %���Ծ���

BD=detJ*transpose(B)*D*B;  %BD����
end

function Kz=StiffnessAssemble(KK,ke,j1,j2,j3,j4,j5,j6,j7,j8)
%ƴװ������󣬸�����Ԫ����ke����ʼ�������KK���ڵ���j��ƴװ�������
%**************��д��������  ������ͨ��ѧ ��еѧԺ ˶1005***************

%��Ԫ������������е�λ��
KLoc=[3*j1-2:3*j1, 3*j2-2:3*j2, 3*j3-2:3*j3, 3*j4-2:3*j4,...
3*j5-2:3*j5, 3*j6-2:3*j6, 3*j7-2:3*j7, 3*j8-2:3*j8];

for ii=1:24 %��Ԫ�ط�����������еĶ�Ӧλ��
for jj=1:24
KK(KLoc(ii),KLoc(jj))=KK(KLoc(ii),KLoc(jj))+ke(ii,jj);
end
end
Kz=KK;
end
