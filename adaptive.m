clear all; 
format long;

%% ========== 区域定义与网格生成 ==========
Lx = 0; Rx = 5;         % x方向边界
Ly = 0; Ry = 1;         % y方向边界
h0 = 0.05;              % 初始网格尺寸

% 定义区域距离函数和尺寸函数
fd = @(p) drectangle(p, Lx, Rx, Ly, Ry);    % 矩形区域距离函数
fh = @(p) h0 + 0.01*p(:,1);                 % 网格尺寸函数（x方向渐变）

% 生成网格
[p, t] = distmesh2d(fd, fh, h0, [Lx, Ly; Rx, Ry], [Lx, Ly; Rx, Ly; Rx, Ry; Lx, Ry]);

%% ========== 邻接关系计算与存储 ==========
fid = fopen('vring_data.m', 'wt');  % 打开邻接数据文件

% 识别边界点索引
indt = [find(p(:,1) < Lx + h0/2)];         % 左边界
indt = [indt; find(p(:,1) > Rx - h0/2)];   % 右边界
indt = [indt; find(p(:,2) < Ly + h0/2)];   % 下边界
indt = [indt; find(p(:,2) > Ry - h0/2)];   % 上边界

id = 1:length(p);       % 所有点初始索引
id(indt) = [];          % 移除边界点（保留内部点）

% 遍历每个内部点构建邻接环
for iter = 1:length(id)
    commonPT = id(iter);                    % 当前中心点
    IndexN = find(t(:) == commonPT);        % 包含该点的三角形索引
    
    % 提取邻接点对
    PT = zeros(length(IndexN), 2);
    for i = 1:length(IndexN)
        % 确定当前三角形顶点
        if mod(IndexN(i), size(t,1)) == 0
            tmpT = t(size(t,1), :);
        else
            tmpT = t(mod(IndexN(i), size(t,1)), :);
        end
        
        % 根据中心点位置确定邻接点对
        if tmpT(1) == commonPT
            PT(i,:) = [tmpT(2), tmpT(3)];
        elseif tmpT(2) == commonPT
            PT(i,:) = [tmpT(3), tmpT(1)];
        else
            PT(i,:) = [tmpT(1), tmpT(2)];
        end
    end
    
    % 对邻接点排序形成闭环
    V = []; 
    V(1:2) = PT(1,:);  % 初始化首两个点
    PT(1,:) = [];      % 移除已处理点
    
    % 按顺序连接邻接点
    for k = 1:length(IndexN)-1
        a = find(PT(:,1) == V(end));
        V(end+1) = PT(a,2);
        PT(a,:) = [];
    end
    
    % 填充至固定长度并写入文件
    V(end+1:8) = 0; 
    fprintf(fid, '%d %d %d %d %d %d %d %d %d \n', length(IndexN), V);
end
fclose(fid);

% 加载邻接数据并验证
D = load('vring_data.m'); 
if max(D(:,1)) > 7
    error('邻接点超过7个，请检查网格');
end

%% ========== 初始条件设置 ==========
valP = zeros(size(p,1), 1);  % 初始化场变量
for i = 1:length(p)
    % 设置初始值（左侧为1，右侧为0）
    if p(i,1) < 0.2
        valP(i) = 1.0;
    else
        valP(i) = 0;
    end
end
nnLap = 0 * valP;  % 拉普拉斯算子初始化

%% ========== 边界条件处理 ==========
% 识别边界点索引
bdid = unique([find(abs(p(:,1)-Lx) < 0.5*h0 | ...
               abs(p(:,1)-Rx) < 0.5*h0 | ...
               abs(p(:,2)-Ly) < 0.5*h0 | ...
               abs(p(:,2)-Ry) < 0.5*h0)]);
           
% 构建边界点邻接关系
nb_in_bd = cell(length(bdid),1);
nb_bd = cell(length(bdid),1);
for i = 1:length(bdid)
    temp = [];
    % 查找共享三角形的所有点
    for j = 1:length(t)
        if sum(bdid(i) == t(j,:)) > 0
            temp = [temp, t(j,:)];
        end
    end
    temp = unique(temp);
    nb_in_bd{i} = [bdid(i), temp];  % 包含边界点的完整环
    
    % 提取内部邻接点
    for k = 1:length(bdid)
        temp(temp == bdid(k)) = [];
    end
    nb_bd{i} = [bdid(i), temp];     % 仅内部邻接点
end

% 设置边界点值为邻接点平均值
for i = 1:length(bdid)
    valP(bdid(i)) = mean(valP(nb_bd{i}(2:end)));
end
nvalP = valP;  % 创建更新副本

%% ========== 空间相关epsilon计算 ==========
deps = valP; 
deps(bdid) = nan;  % 边界点标记为NaN

% 计算每个内部点的局部网格尺寸
epsilon = zeros(size(D,1),1);
for k = 1:size(D,1)
    dd = 0;
    for i = 2:D(k,1)+1
        x0 = p(id(k),1); y0 = p(id(k),2);
        x = p(D(k,i),1); y = p(D(k,i),2);
        dd = dd + sqrt((x0-x)^2 + (y0-y)^2);
    end
    epsilon(k) = 0.5 * dd / D(k,1);  % 平均边长的一半
end
deps(id) = epsilon;  % 更新内部点值

%% ========== 稳定时间步长计算 ==========
dtt = zeros(size(D,1),1);  % 预分配时间步长数组

for k = 1:size(D,1)
    Egrad = 0; 
    Allarea = 0; 
    triC = id(k);  % 中心点
    
    % 遍历邻接三角形
    for i = 2:D(k,1)+1
        % 确定三角形顶点索引
        if i == 2
            triL = D(k, D(k,1)+1);
            triM = D(k,i);
            triN = D(k,i+1);
        else
            triL = D(k,i-1);
            triM = D(k,i);
            triN = D(k,i+1);
        end
        
        % 计算向量和角度
        Ph1 = p(triC,:) - p(triL,:);
        Ph2 = p(triM,:) - p(triL,:);
        Pi1 = p(triC,:) - p(triN,:);
        Pi2 = p(triM,:) - p(triN,:);
        Pj1 = p(triM,:) - p(triC,:);
        Pj2 = p(triN,:) - p(triC,:);
        
        alpA = acos(dot(Ph1, Ph2) / (norm(Ph1)*norm(Ph2)));
        betA = acos(dot(Pi1, Pi2) / (norm(Pi1)*norm(Pi2)));
        
        Egrad = Egrad + (cot(alpA) + cot(betA));  % 梯度权重累加
        Earea = 0.5 * sqrt(norm(Pj1)^2 * norm(Pj2)^2 - dot(Pj1, Pj2)^2); % 三角形面积
        Allarea = Allarea + Earea;  % 总控制面积累加
    end
    
    dtt(k) = 2 * Allarea / (3 * Egrad);  % 局部稳定时间步
end
dt = min(dtt);  % 全局时间步长

%% ========== 主求解循环 ==========
err = 1; 
iter = 0;
while err > 1e-6
    iter = iter + 1;
    nLap = zeros(size(D,1),1);  % 预分配拉普拉斯值
    
    % 遍历每个内部点
    for k = 1:size(D,1)
        Egrad = 0; 
        Allarea = 0; 
        triC = id(k);  % 中心点
        
        % 遍历邻接三角形计算离散拉普拉斯
        for i = 2:D(k,1)+1
            % 确定顶点索引
            if i == 2
                triL = D(k, D(k,1)+1);
                triM = D(k,i);
                triN = D(k,i+1);
            else
                triL = D(k,i-1);
                triM = D(k,i);
                triN = D(k,i+1);
            end
            
            % 计算向量和角度
            Ph1 = p(triC,:) - p(triL,:);
            Ph2 = p(triM,:) - p(triL,:);
            Pi1 = p(triC,:) - p(triN,:);
            Pi2 = p(triM,:) - p(triN,:);
            Pj1 = p(triM,:) - p(triC,:);
            Pj2 = p(triN,:) - p(triC,:);
            
            alpA = acos(dot(Ph1, Ph2) / (norm(Ph1)*norm(Ph2)));
            betA = acos(dot(Pi1, Pi2) / (norm(Pi1)*norm(Pi2)));
            
            % 离散拉普拉斯算子贡献
            Egrad = Egrad + 0.5 * (cot(alpA) + cot(betA)) * (valP(triM) - valP(triC));
            Earea = 0.5 * sqrt(norm(Pj1)^2 * norm(Pj2)^2 - dot(Pj1, Pj2)^2);
            Allarea = Allarea + Earea;  % 控制体积累加
        end
        
        nLap(k) = valP(triC) + dt * 3.0 * Egrad / Allarea;  % 更新拉普拉斯值
    end
    
    % 非线性更新规则
    nvalP(id)=nLap./sqrt(exp(-2*dt./(epsilon.^2))+nLap.^2.*(1-exp(-2*dt./(epsilon.^2))));
    
    % 更新边界值（邻接点平均）
    for bdi = 1:length(bdid)
        nvalP(bdid(bdi)) = mean(nvalP(nb_bd{bdi}(2:end)));
    end
    
    % 计算收敛误差
    err = sqrt(sum((nvalP(id) - valP(id)).^2) / length(id));
    valP = nvalP;  % 更新场变量
    
    % 可视化（每10步）
    if iter == 1 || mod(iter,10) == 0
        fig = figure(1); 
        clf; hold on; 
        view(5,30); 
        box on; grid on;
        fig.Position(3:4) = [1000, 250];
        set(gca, 'FontSize', 13); 
        colormap jet; 
        clim([0,1]);
        
        % 绘制三角形网格
        trimesh(t, p(:,1), p(:,2), nvalP);
        
        % 坐标轴标签
        text([4.5,5.2,-0.3], [-0.05,0.8,0], [-0.15,-0.15,0.75], ...
            {'$x$','$y$','$\phi$'}, 'Interpreter', 'Latex', 'FontSize', 15);
        
        zticks([0,0.5,1]); 
        xticks(0:5); 
        axis([Lx, Rx, Ly, Ry, 0, 1]); 
        drawnow;
    end
end

%% ===== 辅助函数 =====
% 区域距离函数（矩形）
function d = drectangle(p, x1, x2, y1, y2)
    d = -min(min(min(-y1 + p(:,2), y2 - p(:,2)), -x1 + p(:,1)), x2 - p(:,1));
end

% 网格修复函数（省略完整实现以节省空间）
function [p,t,pix]=fixmesh(p,t,ptol)
%FIXMESH  Remove duplicated/unused nodes and fix element orientation.
%   [P,T]=FIXMESH(P,T)

%   Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.

if nargin<3, ptol=1024*eps; end
if nargin>=2 & (isempty(p) | isempty(t)), pix=1:size(p,1); return; end

snap=max(max(p,[],1)-min(p,[],1),[],2)*ptol;
[foo,ix,jx]=unique(round(p/snap)*snap,'rows');
p=p(ix,:);

if nargin>=2
    t=reshape(jx(t),size(t));
    
    [pix,ix1,jx1]=unique(t);
    t=reshape(jx1,size(t));
    p=p(pix,:);
    pix=ix(pix);
    
    if size(t,2)==size(p,2)+1
        flip=simpvol(p,t)<0;
        t(flip,[1,2])=t(flip,[2,1]);
    end
end
end

% 单纯形体积计算（省略完整实现以节省空间）
function v=simpvol(p,t)
%SIMPVOL Simplex volume.
%   V=SIMPVOL(P,T)

%   Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.

switch size(p,2)
 case 1
  d12=p(t(:,2),:)-p(t(:,1),:);
  v=d12;
 case 2
  d12=p(t(:,2),:)-p(t(:,1),:);
  d13=p(t(:,3),:)-p(t(:,1),:);
  v=(d12(:,1).*d13(:,2)-d12(:,2).*d13(:,1))/2;
 case 3
  d12=p(t(:,2),:)-p(t(:,1),:);
  d13=p(t(:,3),:)-p(t(:,1),:);
  d14=p(t(:,4),:)-p(t(:,1),:);
  v=dot(cross(d12,d13,2),d14,2)/6;
 otherwise
  v=zeros(size(t,1),1);
  for ii=1:size(t,1)
    A=zeros(size(p,2)+1);
    A(:,1)=1;
    for jj=1:size(p,2)+1
      A(jj,2:end)=p(t(ii,jj),:);
    end
    v(ii)=det(A);
  end
  v=v/factorial(size(p,2));
end
end

% 网格生成函数（省略完整实现以节省空间）
function [p,t]=distmesh2d(fd,fh,h0,bbox,pfix,varargin)
%DISTMESH2D 2-D Mesh Generator using Distance Functions.
%   [P,T]=DISTMESH2D(FD,FH,H0,BBOX,PFIX,FPARAMS)
%
%      P:         Node positions (Nx2)
%      T:         Triangle indices (NTx3)
%      FD:        Distance function d(x,y)
%      FH:        Scaled edge length function h(x,y)
%      H0:        Initial edge length
%      BBOX:      Bounding box [xmin,ymin; xmax,ymax]
%      PFIX:      Fixed node positions (NFIXx2)
%      FPARAMS:   Additional parameters passed to FD and FH
%
%   Example: (Uniform Mesh on Unit Circle)
%      fd=@(p) sqrt(sum(p.^2,2))-1;
%      [p,t]=distmesh2d(fd,@huniform,0.2,[-1,-1;1,1],[]);
%
%   Example: (Rectangle with circular hole, refined at circle boundary)
%      fd=@(p) ddiff(drectangle(p,-1,1,-1,1),dcircle(p,0,0,0.5));
%      fh=@(p) 0.05+0.3*dcircle(p,0,0,0.5);
%      [p,t]=distmesh2d(fd,fh,0.05,[-1,-1;1,1],[-1,-1;-1,1;1,-1;1,1]);
%
%   Example: (Polygon)
%      pv=[-0.4 -0.5;0.4 -0.2;0.4 -0.7;1.5 -0.4;0.9 0.1;
%          1.6 0.8;0.5 0.5;0.2 1;0.1 0.4;-0.7 0.7;-0.4 -0.5];
%      [p,t]=distmesh2d(@dpoly,@huniform,0.1,[-1,-1; 2,1],pv,pv);
%
%   Example: (Ellipse)
%      fd=@(p) p(:,1).^2/2^2+p(:,2).^2/1^2-1;
%      [p,t]=distmesh2d(fd,@huniform,0.2,[-2,-1;2,1],[]);
%
%   Example: (Square, with size function point and line sources)
%      fd=@(p) drectangle(p,0,1,0,1);
%      fh=@(p) min(min(0.01+0.3*abs(dcircle(p,0,0,0)), ...
%                   0.025+0.3*abs(dpoly(p,[0.3,0.7; 0.7,0.5]))),0.15);
%      [p,t]=distmesh2d(fd,fh,0.01,[0,0;1,1],[0,0;1,0;0,1;1,1]);
%
%   Example: (NACA0012 airfoil)
%      hlead=0.01; htrail=0.04; hmax=2; circx=2; circr=4;
%      a=.12/.2*[0.2969,-0.1260,-0.3516,0.2843,-0.1036];
%
%      fd=@(p) ddiff(dcircle(p,circx,0,circr),(abs(p(:,2))-polyval([a(5:-1:2),0],p(:,1))).^2-a(1)^2*p(:,1));
%      fh=@(p) min(min(hlead+0.3*dcircle(p,0,0,0),htrail+0.3*dcircle(p,1,0,0)),hmax);
%
%      fixx=1-htrail*cumsum(1.3.^(0:4)');
%      fixy=a(1)*sqrt(fixx)+polyval([a(5:-1:2),0],fixx);
%      fix=[[circx+[-1,1,0,0]*circr; 0,0,circr*[-1,1]]'; 0,0; 1,0; fixx,fixy; fixx,-fixy];
%      box=[circx-circr,-circr; circx+circr,circr];
%      h0=min([hlead,htrail,hmax]);
%
%      [p,t]=distmesh2d(fd,fh,h0,box,fix);

%
%   See also: MESHDEMO2D, DISTMESHND, DELAUNAYN, TRIMESH.

%   distmesh2d.m v1.1
%   Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.

dptol=.001; ttol=.1; Fscale=1.2; deltat=.2; geps=.001*h0; deps=sqrt(eps)*h0;
densityctrlfreq=30;

% 1. Create initial distribution in bounding box (equilateral triangles)
[x,y]=meshgrid(bbox(1,1):h0:bbox(2,1),bbox(1,2):h0*sqrt(3)/2:bbox(2,2));
x(2:2:end,:)=x(2:2:end,:)+h0/2;                      % Shift even rows
p=[x(:),y(:)];                                       % List of node coordinates

% 2. Remove points outside the region, apply the rejection method
p=p(feval(fd,p,varargin{:})<geps,:);                 % Keep only d<0 points
r0=1./feval(fh,p,varargin{:}).^2;                    % Probability to keep point
p=p(rand(size(p,1),1)<r0./max(r0),:);                % Rejection method
if ~isempty(pfix), p=setdiff(p,pfix,'rows'); end     % Remove duplicated nodes
pfix=unique(pfix,'rows'); nfix=size(pfix,1);
p=[pfix; p];                                         % Prepend fix points
N=size(p,1);                                         % Number of points N

count=0;
pold=inf;                                            % For first iteration
% clf,view(2),axis equal,axis off
while 1
  count=count+1;
  % 3. Retriangulation by the Delaunay algorithm
  if max(sqrt(sum((p-pold).^2,2))/h0)>ttol           % Any large movement?
    pold=p;                                          % Save current positions
    t=delaunayn(p);                                  % List of triangles
    pmid=(p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;    % Compute centroids
    t=t(feval(fd,pmid,varargin{:})<-geps,:);         % Keep interior triangles
    % 4. Describe each bar by a unique pair of nodes
    bars=[t(:,[1,2]);t(:,[1,3]);t(:,[2,3])];         % Interior bars duplicated
    bars=unique(sort(bars,2),'rows');                % Bars as node pairs
    % 5. Graphical output of the current mesh
%     cla,patch('vertices',p,'faces',t,'edgecol','k','facecol',[.8,.9,1]);
%     drawnow
  end

  % 6. Move mesh points based on bar lengths L and forces F
  barvec=p(bars(:,1),:)-p(bars(:,2),:);              % List of bar vectors
  L=sqrt(sum(barvec.^2,2));                          % L = Bar lengths
  hbars=feval(fh,(p(bars(:,1),:)+p(bars(:,2),:))/2,varargin{:});
  L0=hbars*Fscale*sqrt(sum(L.^2)/sum(hbars.^2));     % L0 = Desired lengths
  
  % Density control - remove points that are too close
  if mod(count,densityctrlfreq)==0 & any(L0>2*L)
      p(setdiff(reshape(bars(L0>2*L,:),[],1),1:nfix),:)=[];
      N=size(p,1); pold=inf;
      continue;
  end
  
  F=max(L0-L,0);                                     % Bar forces (scalars)
  Fvec=F./L*[1,1].*barvec;                           % Bar forces (x,y components)
  Ftot=full(sparse(bars(:,[1,1,2,2]),ones(size(F))*[1,2,1,2],[Fvec,-Fvec],N,2));
  Ftot(1:size(pfix,1),:)=0;                          % Force = 0 at fixed points
  p=p+deltat*Ftot;                                   % Update node positions

  % 7. Bring outside points back to the boundary
  d=feval(fd,p,varargin{:}); ix=d>0;                 % Find points outside (d>0)
  dgradx=(feval(fd,[p(ix,1)+deps,p(ix,2)],varargin{:})-d(ix))/deps; % Numerical
  dgrady=(feval(fd,[p(ix,1),p(ix,2)+deps],varargin{:})-d(ix))/deps; %    gradient
  dgrad2=dgradx.^2+dgrady.^2;
  p(ix,:)=p(ix,:)-[d(ix).*dgradx./dgrad2,d(ix).*dgrady./dgrad2];    % Project

  % 8. Termination criterion: All interior nodes move less than dptol (scaled)
  if max(sqrt(sum(deltat*Ftot(d<-geps,:).^2,2))/h0)<dptol, break; end
end

% Clean up and plot final mesh
[p,t]=fixmesh(p,t);
% simpplot(p,t)
end