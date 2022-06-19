function [solveVec,solveMatrix]=Model17()
%MODEL17 此处显示有关此函数的摘要
%   此处显示详细说明

%%%%%%%%%%%%%%%%%%%%I参数更改区%%%%%%%%%%%%%%%%%%%%%%%%
k=100;       %直角边划分几个点
spanT=800;  %非稳态计算(spanT-1)次,
            % 即[0,delta_t,...,(spanT-1)delta_t]共spanT-1张图片
a=0.3;  %三角形x方向[0,a]
b=0.3;  %三角形y方向[0,b]
delta_t=0.005;  %计算非稳态传热的间隔时间
filename='temperature01.gif';   %保存的文件名
delay_time=0.2; %gif各帧的延迟时间
Vec1=[1,160,0,0]; %边1为第2类边值条件,热流强度300
% Vec1=[1,300,0,0]; %边1为第2类边值条件,热流强度300
% Vec1=[2,0,20,100];  %边1为第3类边值条件,表面传热系数20,外界环境温度100
Vec2=[2,0,20,100];  %边2为第3类边值条件,表面传热系数20,外界环境温度100
Vec3=[2,0,20,20]; %边3为第3类边值条件,表面传热系数20,外界环境温度20
% Vec3=[2,0,20,100];  %边3为第3类边值条件,表面传热系数20,外界环境温度20
Vec4=[0,0,0,0];     %内部状态
ka=1;   %热传导系数
lou=5;  %密度
ca=20;  %比热容
Pos=[ 0.06 0.11 0.37 0.76;
    0.52 0.11 0.42 0.76 ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%II建模与数据处理%%%%%%%%%%%%%%%%%%%%%
disp('data processing...');
delta_x=a/(k-1);
delta_y=b/(k-1);
loc_data=zeros(2,k*(k+1)/2);

% lou=5;
% ca=20;

rows=@(i,k) (2*k+2-i)*(i-1)/2;
side1=zeros(k,1);
side2=zeros(k,1);
side3=zeros(k,1);

count=1;
for i=1:k
    loc_y=(i-1)*delta_y;
    for j=1:k+1-i
        if i==1
            side2(count)=count;
        end
        if j==1
            side1(count)=count;
        end
        if j==k+1-i
            side3(count)=count;
        end
        loc_x=(j-1)*delta_x;
        loc_data(:,count)=[loc_x;loc_y];
        count=count+1;
    end
end

side1=side1(side1~=0);
side2=side2(side2~=0);
side3=side3(side3~=0);


numMatrix=zeros(3,(k-1)^2);
count=1;
for i=2:k
    for j=1:k-i
        numMatrix(:,count)=[rows(i,k)+j,rows(i-1,k)+j,rows(i-1,k)+j+1]';
        numMatrix(:,count+1)=[rows(i,k)+j,rows(i,k)+j+1,rows(i-1,k)+j+1]';
        count=count+2;
    end
    j=k+1-i;
    numMatrix(:,count)=[rows(i,k)+j,rows(i-1,k)+j,rows(i-1,k)+j+1]';
    count=count+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%III矩阵构建于求解%%%%%%%%%%%%%%%%%%%%%%%%
% ka=1;
Kmatrix=zeros(size(loc_data,2),size(loc_data,2));
Cmatrix=zeros(size(loc_data,2),size(loc_data,2));
Kvec=zeros(size(loc_data,2),1);
% 在后面再转成sparse
% Kmatrix=sparse(size(loc_data,2),size(loc_data,2));
% Cmatrix=sparse(size(loc_data,2),size(loc_data,2));

% Vec1=[1,160,0,0];
% Vec2=[2,0,20,100];
% Vec3=[2,0,20,20];
% Vec4=[0,0,0,0];
infVec=cell(3,1);
indexer=1:3;

for i=1:size(numMatrix,2)
    infVec(1)={Vec4};
    infVec(2)={Vec4};
    infVec(3)={Vec4};
    serial=numMatrix(:,i);

    judge1=[any(side1==serial(1)),any(side1==serial(2)),any(side1==serial(3))];
    if sum(judge1)==2
        infVec(mod(indexer(judge1==false),3)+1)={Vec1};
    end

    judge2=[any(side2==serial(1)),any(side2==serial(2)),any(side2==serial(3))];
    if sum(judge2)==2
        infVec(mod(indexer(judge2==false),3)+1)={Vec2};
    end

    judge3=[any(side3==serial(1)),any(side3==serial(2)),any(side3==serial(3))];
    if sum(judge3)==2
        infVec(mod(indexer(judge3==false),3)+1)={Vec3};
    end
    [Kmatrix,Kvec,Cmatrix]=subModel(serial,infVec{1},infVec{2},infVec{3},loc_data,ka,Kmatrix,Kvec,lou,ca,Cmatrix);
end



disp('matrix constructed')
disp('solving linear systems...');

% spanT=800;
solveMatrix=zeros(size(loc_data,2),spanT);
solveVec=Kmatrix\Kvec;

% delta_t=0.005;
solveMatrix(:,1)=50;
Cmatrix=sparse(Cmatrix);
Kmatrix=sparse(Kmatrix);
Smatrix=Cmatrix+0.5*delta_t*Kmatrix;
for m=2:spanT
    solveMatrix(:,m)=Smatrix\(Cmatrix*solveMatrix(:,m-1)+0.5*delta_t*(2*Kvec-Kmatrix*solveMatrix(:,m-1)));
    % solveMatrix(:,m)=(Cmatrix+0.5*delta_t*Kmatrix)\(Cmatrix*solveMatrix(:,m-1)+0.5*delta_t*(2*Kvec-Kmatrix*solveMatrix(:,m-1)));
end

disp('solved')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%IV可视化%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('visualizing...')

rangeMax=max([solveVec,solveMatrix],[],'all')+1e-6;
rangeMin=min([solveVec,solveMatrix],[],'all')-1e-6;
range=linspace(rangeMin,rangeMax,21);

close all;
set(gcf,'visible','on','Position',[15,100,1400,650]);

pic_num=1;

Xmatrix=zeros(3,(k-1)^2);
Ymatrix=zeros(3,(k-1)^2);
Tmatrix=zeros(3,(k-1)^2);

colormap(jet(length(range)-1));
c=colorbar;
c.TickLabels=range;
c.Ticks=linspace(0,1,length(range));
axis equal;


for i=1:(k-1)^2
    numarr=numMatrix(:,i);    
    Xmatrix(:,i)=loc_data(1,numarr)';
    Ymatrix(:,i)=loc_data(2,numarr)';
    Tmatrix(:,i)=solveVec(numarr);
end

ColorsmZero=colorMaper2(mean(Tmatrix),range);


for m=1:spanT

    hy=subplot(1,2,1);
    hy.Position=Pos(1,:);

    % colorbar('Ticks',[0:1/(length(range)-1):1],'TickLabels',range);
    title('稳态温度场');
    % colorbar('Ticks',[0:1/(length(range)-1):1],'TickLabels',range);
    patch(Xmatrix,Ymatrix,ColorsmZero);
    Tmatrix=reshape(solveMatrix(numMatrix,m),[3,(k-1)^2]);
    Colorsm=colorMaper2(mean(Tmatrix),range);
    hold on;

    hy=subplot(1,2,2);
    hy.Position=Pos(2,:);

    tt=num2str(mod((m-1)*delta_t,1));
    if mod(m,2)==1 && length(tt)>2
        tail=repmat('0',[1,5-length(tt)]);
        title(['温度场',num2str(m),':',num2str((m-1)*delta_t),tail,'s']);
    elseif length(tt)<2
        tail=repmat('0',[1,3]);
        title(['温度场',num2str(m),':',num2str((m-1)*delta_t),'.',tail,'s']);
    else
        title(['温度场',num2str(m),':',num2str((m-1)*delta_t),'s']);
    end

    patch(Xmatrix,Ymatrix,Colorsm);
    c=colorbar;
    c.TickLabels=range;
    c.Ticks=linspace(0,1,length(range));

    disp(['figure:',num2str(m)]);
    hold on;
    drawnow
    F=getframe(gcf);
    I=frame2im(F);
    [I,map]=rgb2ind(I,256);
    if pic_num == 1
        imwrite(I,map,filename,'gif','Loopcount',inf,'DelayTime',delay_time);
    else
        imwrite(I,map,filename,'gif','WriteMode','append','DelayTime',delay_time);
    end
    pic_num = pic_num + 1;
    if m<spanT
        clf
    end

end
disp('finished');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

function result=colorMaper2(Tvec,range)
n=length(Tvec);
step=length(range)-1;
color=jet(step);
indexer=(1:step)';
tempM=range'-Tvec;
tt=sum(indexer.*(tempM(1:step,:).*tempM(2:step+1,:)<=0));
if any(tt==0)
    indexError=find(tt==0);
    sign=tempM(2,tt==0)>0;
    tt(indexError(sign))=step;
    tt(indexError(~sign))=1;
end
if any(tt>step)
    indexError=find(tt==0);
    tt(indexError)=floor(tt(indexError)/2);
end
result=reshape(color(tt,:),[n,1,3]);
end

function [Kmatrix,Kvec,Cmatrix]=subModel(serial,infVec1,infVec2,infVec3,loc,ka,Kmatrix,Kvec,lou,ca,Cmatrix)
%SUBMODEL 此处显示有关此函数的摘要
%   此处显示详细说明
%   serial=[number1,number2,number3]
%   infVec是4*1s的信息块
%   infVec(1):类型
%   infVec(2):q_f
%   infVec(3):h_c
%   infVec(4):T_0

num1=serial(1);
num2=serial(2);
num3=serial(3);
tempM=[2,1;1,2];
tempM2=[2,1,1;1,2,1;1,1,2];

space=abs(0.5*det([loc(:,num2),loc(:,num3)]-loc(:,num1)));

Minv=inv([loc(:,num1)',1;loc(:,num2)',1;loc(:,num3)',1]);

Kmatrix([num1,num2,num3],[num1,num2,num3])=Kmatrix([num1,num2,num3],[num1,num2,num3])+ka*space*Minv(1:2,:)'*Minv(1:2,:);

%num1-num2
if infVec1(1)==0
elseif infVec1(1)==1
    side=norm(loc(:,num2)-loc(:,num1));
    Kvec([num1,num2])=Kvec([num1,num2])+infVec1(2)*side*0.5;
elseif infVec1(1)==2
    side=norm(loc(:,num2)-loc(:,num1));
    Kmatrix([num1,num2],[num1,num2])=Kmatrix([num1,num2],[num1,num2])+infVec1(3)*side/6*tempM;
    Kvec([num1,num2])=Kvec([num1,num2])+infVec1(3)*infVec1(4)*side*0.5;
end

%num2-num3
if infVec2(1)==0
elseif infVec2(1)==1
    side=norm(loc(:,num3)-loc(:,num2));
    Kvec([num2,num3])=Kvec([num2,num3])+infVec2(2)*side*0.5;
elseif infVec2(1)==2
    side=norm(loc(:,num3)-loc(:,num2));
    Kmatrix([num2,num3],[num2,num3])=Kmatrix([num2,num3],[num2,num3])+infVec2(3)*side/6*tempM;
    Kvec([num2,num3])=Kvec([num2,num3])+infVec2(3)*infVec2(4)*side*0.5;
end

%num3-num1
if infVec3(1)==0
elseif infVec3(1)==1
    side=norm(loc(:,num3)-loc(:,num1));
    Kvec([num3,num1])=Kvec([num3,num1])+infVec3(2)*side*0.5;
elseif infVec3(1)==2
    side=norm(loc(:,num3)-loc(:,num1));
    Kmatrix([num3,num1],[num3,num1])=Kmatrix([num3,num1],[num3,num1])+infVec3(3)*side/6*tempM;
    Kvec([num1,num3])=Kvec([num1,num3])+infVec3(3)*infVec3(4)*side*0.5;
end

Cmatrix([num1,num2,num3],[num1,num2,num3])=Cmatrix([num1,num2,num3],[num1,num2,num3])+lou*ca*space/12*tempM2;



end







