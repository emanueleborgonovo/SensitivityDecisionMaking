function [y,dltas,flatcnt]=partialdep2(x,f,varargin) %$ ,lb,ub,lbls
% PARTIALDEP Partial Dependence Plots.
% Y=PARTIALDEP(X,F) returns the partial dependence curves of F(X),
%    producing a graphical output. If F is a vector, a regression model is
%    invoked.
% Y=PARTIALDEP(X,F,'OPT',ARG,...) offers lots of undocumented features.
[n,k]=size(x);
% default options
opts=struct('GfxTitle','',...
            'GfxCols',min(k,3),...
            'Labels','',...
            'OutputLabel','y',...
            'Inputs',1:k,...
            'DiscreteInputs',[],...
	        'LowerBound',min(x),...
            'UpperBound',max(x),...
            'Points',100,...
            'ShowEnvelop',false,...
            'ShowQuantiles',false,...
            'ShowOneWay',30,...
            'ShowScatter',false,...
            'GradientInfo',[],...
            'ScatterData',[],...
            'SecondOrder',[],...
			'SpiderPlot',false,...
            'InputConstraints',@(x)true(size(x,1),1),...
            'Light',[]); %,...
%            'ResponseModel',@fitrtree);
%
%disp('Using the local test version of partialdep.');
%
if(nargin>=3) && ~isempty(varargin)
    if isstruct(varargin{1})
        opts_in=varargin{1};
    else
        opts_in=struct(varargin{:});
    end
    members=fieldnames(opts);
    for i=1:length(members)
            o=members{i};
            if isfield(opts_in,o)
                opts.(o)=opts_in.(o);
            end
    end
end
l=length(opts.Inputs);
if(opts.SpiderPlot)
    if(~isempty(opts.DiscreteInputs))
        warning('No spider plot for discrete inputs.');
        opts.Spiederplot=false;
    end
end
if(~isempty(opts.GradientInfo))
 js=sort(randperm(opts.Points,opts.GradientInfo));
 dltas=nan(1,k);
 flatcnt=zeros(1,k);dflat=2/(opts.GradientInfo*(opts.GradientInfo-1));
else 
 js=0;
end

if(isa(f,'function_handle') || contains(class(f),'Regression'))
    %% regression model provided
m=1;
    for i=opts.Inputs
    yys=[]; % for GradientInfo
      if(~ismember(i,opts.DiscreteInputs))
       OneWay=min(opts.ShowOneWay,opts.Points);
       ii=randperm(opts.Points,OneWay);
       z=linspace(opts.LowerBound(i),opts.UpperBound(i),opts.Points);
       w=nan(opts.Points,opts.ShowOneWay);
       plt_str='-';
      else
          % for discrete, only use available data
          z=unique(x(:,i))';
          OneWay=min(opts.ShowOneWay,length(z));
          ii=randperm(opts.Points,OneWay);
          w=nan(opts.Points,length(z));
          plt_str='-o';
      end
      % sizes may change between iterations (with discrete inputs)
      clear p q r q4 w
      for j=1:length(z)
        xx=x;xx(:,i)=z(j); %*ones(n,1);
        support=opts.InputConstraints(xx);
        allsupport=all(support);
        jj=cumsum(support);
        xx(~support,:)=[];
        if(isa(f,'function_handle'))
            yy=f(xx);
        else
            yy=predict(f,xx);
        end
        p(j)=mean(yy);
        % save for Gradient Info, may fail when constraints are active
        if any(j==js), yys=[yys,yy]; end
        %if(j==j1), yy1=yy; end
        %if(opts.ShowEnvelop)
        q(j)=min([yy;inf]);
        r(j)=max([yy;-inf]);
        %end
        if(opts.ShowQuantiles)
            if(islogical(opts.ShowQuantiles)), Q=4; else
                Q=opts.ShowQuantiles; end
            ys=sort([-inf;yy;inf]);
            q4(j,:)=ys(ceil((n+2)*(1:(Q-1))/Q));
        end
        if(opts.ShowOneWay)
            if(allsupport)
                w(j,:)=yy(ii);
            else
               % alignment due to removed data 
                for mm=1:length(ii)
                    if(support(ii(mm))) 
                         w(j,mm)=yy(jj(ii(mm))); 
                    end
                end
            end   
        end
    end
	if(~opts.SpiderPlot)
     subplot(ceil(l/opts.GfxCols),opts.GfxCols,m);m=m+1;
     if(opts.ShowScatter && ~isempty(opts.ScatterData))
        plot(x(:,i),opts.ScatterData,'b*','MarkerSize',1);
        hold on;
     end
     if(opts.ShowOneWay)
        plot(z,w,['r' plt_str],'LineWidth',1);
        hold on;
     end
     plot(z,p,['k' plt_str],'LineWidth',3);
%    a=axis;
     a(1)=opts.LowerBound(i);
     a(2)=opts.UpperBound(i);
     a(3)=min(q);a(4)=max(r);
     axis(a);
     if(~isempty(opts.GradientInfo))
      if(opts.Points==length(z)) % 
        thresh0=(max(yy)-min(yy))*0.01; % constant mode 
        dlta=zeros(opts.GradientInfo*(opts.GradientInfo-1)/2,1);
        kk=1;
        for ll=1:(opts.GradientInfo-1)
            for mm=(ll+1):opts.GradientInfo
            g0=p(js(mm))-p(js(ll));
            g1=yys(:,mm)-yys(:,ll);
            
            if abs(g0)>thresh0
                dlta(kk)=mean(sign(g1./g0));
            else
                dlta(kk)=mean(2*(abs(g1)<thresh0)-1); % g1==0
                flatcnt(i)=flatcnt(i)+dflat;
            end
            kk=kk+1;
            end
        end
        dltas(i)=.5-.5*mean(dlta); 
      else
        flatcnt(i)=-99; % magic number    
      end
     end
     if(opts.ShowEnvelop)
        hold on;
        plot(z,q,'r:',z,r,'r:','LineWidth',2);
     end
     if(opts.ShowQuantiles)
        hold on;
        plot(z,q4,'b:','LineWidth',1);
     end
	 hold off
	 if~isempty(opts.GfxTitle), title(opts.GfxTitle);
     elseif(~isempty(opts.GradientInfo))
      if(flatcnt(i)~=-99)
       title(['Discrepancy ',num2str(dltas(i),'%0.2f'),...
          ', Flatness ',num2str(flatcnt(i),'%0.2f')])
      else
          title('No Gradient Info')
      end
     end 
	 if(isempty(opts.Labels))
      xlabel(['x_{' num2str(i) '}']);
	 else
      xlabel(opts.Labels{i},'Interpreter','none');
	 end
     ylabel(opts.OutputLabel);
     set(gca,'FontSize',14);
    end
    if(opts.SpiderPlot)
     y(i,:)=p;
    end
end
if(opts.SpiderPlot)
 if(any(opts.UpperBound(opts.Inputs).*opts.LowerBound(opts.Inputs)<0)), error('Sign change detected.'); end

 xs=opts.LowerBound(opts.Inputs)+ bsxfun(@times,opts.UpperBound(opts.Inputs)-opts.LowerBound(opts.Inputs),linspace(0,1,100)');
 plot(bsxfun(@rdivide,xs,mean(x(:,opts.Inputs))),y');
 xlabel('Change from mean case');
     ylabel(opts.OutputLabel);
 %set(gca,'FontSize',14);
 if(~isempty(opts.Labels))
     h=legend(opts.Labels,'AutoUpdate','off','NumColumns',2);
	 end
   title(opts.GfxTitle);
% Percentages on the X-Axis
     a=[cellstr(num2str(round(get(gca,'xtick')'*100-100),'%+d'))]; 
     pct = char(ones(size(a,1),1)*'%'); 
     set(gca,'xticklabel',[char(a),pct]);
end
else
 if(k>5)
    % use regression tree implementation
    mdl=fitrtree(x,f);
else
    % use gaussian process regression
     mdl=fitrgp(x,f,'Standardize',true,'Verbose',1);
    % use support vector machine regression
    %mdl=fitrsvm(x,f);
 end
 % call self with regression model argument
 y=partialdep(x,mdl,'ScatterData',f,varargin{:});
%  m=1;
%     for i=opts.Inputs
%         subplot(ceil(l/opts.GfxCols),opts.GfxCols,m);m=m+1;
%         if(opts.ShowScatter)
%             plot(x(:,i),y,'b*','MarkSize',2);
%             hold on;
%         end
%         plotPartialDependence(mdl,i);
%     end
end
if(~isempty(opts.SecondOrder))
 figure
 i=opts.SecondOrder(1);j=opts.SecondOrder(2);
 p2=ceil(sqrt(opts.Points));
 [xxi,xxj]=meshgrid(linspace(opts.LowerBound(i),opts.UpperBound(i),p2),...
    linspace(opts.LowerBound(j),opts.UpperBound(j),p2));
 
 for l=1:numel(xxi)
    xx=x;
    xx(:,i)=xxi(l);
    xx(:,j)=xxj(l);
    xx(~opts.InputConstraints(xx),:)=[];

    if(exist('mdl'))
        yy=predict(mdl,xx);
    else
    if(isa(f,'function_handle'))
        yy=f(xx);
    else
        yy=predict(f,xx);
    end;end
    pp(l)=mean(yy);
 end
 
% mesh(xxi,xxj,reshape(pp',p2,p2));
surf(xxi,xxj,reshape(pp',p2,p2),'FaceColor','interp','FaceLighting','gouraud');
if~isempty(opts.Light), light('Position',opts.Light);end
if(isempty(opts.Labels))
    xlabel(['x_{' num2str(i) '}']);
    ylabel(['x_{' num2str(j) '}']);
else
    xlabel(opts.Labels{i});
    ylabel(opts.Labels{j});
end
end
end

function testpartialdep
%%
    n0=2^15;k=2;
u=sobolseq(n0,k);x=u(sum(u.^2,2)>.5,:); % dependence
n=size(x,1);
model=@(x)x(:,1).^2+x(:,2)/6;
y=model(x);
%%
partialdep(x,model);
%%
    GfxRows=3;
for i=1:k,subplot(ceil(k/GfxRows),GfxRows,i);hold on; end
%% pause
partialdep(x,y);
%%
%p=@(x)(sum(x.^2,2)>.5) ./(1-pi/8);
%pj=@(x)((1-(x<sqrt(.5)).*sqrt(.5-x.^2)))./(1-pi/8);
%
%cosi(x,y.*pj(x(:,1)).*pj(x(:,2))./p(x),8,'weighted');
end

function testwithconstraints
%%
clf
mdl=@(x)x(:,1)+x(:,2);
x=2*rand(1000,2)-1;
constr=@(x)(x(:,1).^2+x(:,2).^2)<=1;
x(~constr(x),:)=[];
partialdep(x,mdl,'ShowOneWay',30);
pause(5);clf
partialdep(x,mdl,'InputConstraints',constr,'ShowOneWay',30)
pause(5);
clf
x=2*rand(10000,2)-1;
constr=@(x)(x(:,1)+x(:,2))<=0;
x(~constr(x),:)=[];
partialdep(x,mdl,'ShowOneWay',30);
pause(5);clf
partialdep(x,mdl,'InputConstraints',constr,'ShowOneWay',30)
pause(5);clf
x=2*rand(10000,2)-1;
constr=@(x)((2*x(:,1)+x(:,2))<=-1) | ((x(:,1)+2*x(:,2))>=1)
x(~constr(x),:)=[];
partialdep(x,mdl,'ShowOneWay',30);
pause(5);clf
partialdep(x,mdl,'InputConstraints',constr,'ShowOneWay',30)
pause(5);clf;
x=rand(10000,2);
mdl=@(x)sum(x,2);
constr=@(x)((x(:,2))<=1-2*x(:,1)) | ((x(:,2))>=1-x(:,1)/2);
x(~constr(x),:)=[];
partialdep(x,mdl,'InputConstraints',constr,'ShowOneWay',30);
%%
end