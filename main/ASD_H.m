%Amplitude Spetral Density Helper version 1.0.7 (Marmota)
%
%Author email: c.houyuan@mail.scut.edu.cn
%
%用于分析计算离散时序信号X(t_n)的均方根功率谱密度ASD，用法如下：
%
% -- [freq, ASD] = ASD_H(X, dt);
% -- [freq, ASD] = ASD_H(X, dt, Win_func);
% -- [freq, ASD] = ASD_H(X, dt, Win_func, alpha);
% -- [freq, ASD] = ASD_H(X, dt, Win_Length);
% -- [freq, ASD] = ASD_H(X, dt, Win_Length, overlap);
% -- [freq, ASD] = ASD_H(X, dt, Win_Length, overlap, Win_func);
% -- [freq, ASD] = ASD_H(X, dt, Win_Length, overlap, Win_func, plot_option);
% -- [freq, ASD] = ASD_H(X, dt, Win_Length, overlap, Win_func, alpha, plot_option);
%
%
%    X, dt	一维数组及其时间间隔，这两个参数作为时序信号的基本输入。
%
%    Win_Length	窗函数的时域宽度，不能设置超过原信号X的时域宽度，以秒为单位。
%		默认情况下窗宽会自动设置为原信号的总时长。
%		比如原信号总长度为1小时，那么默认Win_Length=3600。
%
%    overlap	窗函数的重叠度，默认为不重叠，即overlap=0，
%		注意该参数只能取0.0-0.9之间的值。
%
%    Win_func	指定窗函数的类型，默认选择Win_func="hanning"，目前支持的窗函数：
%
%	"hanning", "hann", "hn" (default)
%		ENBW			 1.50
%		Sidelobe Falloff	-20.7 dB/octave, -68.9 dB/decade
%
%	"hann-poisson", "hp"		Hann-Poisson Window
%		ENBW			 1.61
%		Sidelobe Falloff	-16.8 dB/octave, -55.9 dB/decade
%
%	"kaiser", "ks" (alpha=10)	Kaiser Window
%		ENBW			 1.76
%		Sidelobe Falloff	-50+ dB/octave, -100+ dB/decade
%
%	"kaiser-bessel", "kb"		Kaiser-Bessel Window
%		ENBW			 1.80
%		Sidelobe Falloff	-11.8 dB/octave, -39.3 dB/decade
%
%	"gauss", "norm" (alpha=3.33)	Gaussian Window
%		ENBW			 1.89
%		Sidelobe Falloff	-10.5 dB/octave, -34.8 dB/decade
%
%	"cgauss", "cnorm" (alpha=2)	Confined Gaussian Window
%		ENBW			 1.32
%		Sidelobe Falloff	-16.3 dB/octave, -54.3 dB/decade
%
%	"pcos" (alpha=3)		Power of Cosine Window
%		ENBW			 1.74
%		Sidelobe Falloff	-21.3 dB/octave, -70.6 dB/decade
%
%	"hamming", "hamm", "hm"
%		ENBW			 1.36
%		Sidelobe Falloff	-7.90 dB/octave, -26.1 dB/decade		
%
%	"flattop", "ft"
%		ENBW			 3.78
%		Sidelobe Falloff	-12.1 dB/octave, -40.1 dB/decade
%
%	"blackman", "bm"
%		ENBW			 1.73
%		Sidelobe Falloff	-21.6 dB/octave, -71.7 dB/decade
%
%	"bohman"
%		ENBW			 1.79
%		Sidelobe Falloff	-21.4 dB/octave, -71.0 dB/decade
%
%	"triangular", "bartlett", "tr"	Triangle Window
%		ENBW			 1.34
%		Sidelobe Falloff	-11.4 dB/octave, -37.7 dB/decade
%
%	"exp", "exponential", "poisson"
%		ENBW			 1.07
%		Sidelobe Falloff	-6.50 dB/octave, -21.8 dB/decade
%
%	"welch"
%		ENBW			 1.20
%		Sidelobe Falloff	-11.0 dB/octave, -36.5 dB/decade
%
%	"nuttall","nt"
%		ENBW			 2.03
%		Sidelobe Falloff	-23.4 dB/octave, -77.6 dB/decade
%
%	"planck-taper","pt" (alpha=2.5)
%		ENBW			 2.03
%		Sidelobe Falloff	-30.9 dB/octave, -102.6 dB/decade
%
%	"rectangular", "rect", "none", "box" (not recommended)
%		ENBW			 1.0
%		Sidelobe Falloff	-6.00 dB/octave, -20.0 dB/decade
%
%	Template reference:
%	https://www.recordingblogs.com/wiki/window
%
%
%    plot_option:
%	"r", "b", "g"	分别代表红蓝绿三种线条颜色，
%	"-", "-.", "--", ":."	线型
%	".", "*", "o", "^"	点图
%	实际绘图规则参见函数plot()，使用范例：
%
%	  o均方根功率谱和频率矢量输出的同时绘制红色虚线预览图
%	  -- [freq, ASD] = ASD_H (..., "r--");
%
%	  o仅绘制蓝色点阵预览图
%	  -- ASD_H (..., "b.");
%
%
%   输出：
%   freq	时域下的频率矢量，以Hz为量纲；
%    ASD	均方根功率谱密度矢量，与freq逐一对应，以X/sqrt(Hz)为量纲，这里X意思是时序数据X(t)本身的量纲。
%
%Marmota本本更新:
%    新增平滑模式，使用示例如下，
%
% -- 假设时序数据为X，采样时间为60s，这里推荐使用kaiser窗（alpha参数默认为10），
%
%	使用默认的平滑参数，保留低频的个别点，在高频按对数分区并求平均
% -- [freq, ASD] = ASD_H(X, 60, "ks", 10, "smooth", "b-");
%
%	使用默认的平滑参数，保留低频的200个点，在高频按对数作1000个分区并求平均
% -- [freq, ASD] = ASD_H(X, 60, "ks", 10, "smooth", [200 1000], "b-");
%


function [freq,ASD]=ASD_H(X, dt, varargin)

Win_Spec="hanning";
L=length(X);
firstnumeric=1;quickview=0;Win_recov=2;alpha_specified=0;smooth=0;smpara_specified=0;
Twin_specified=0;overlap_specified=0;smpara=[];skip=0;nidx=[];
T_winWidth=(L-1)*dt;
overlap=0;
tn=(1:L);
if size(X,1)>1
tn=tn';
end
coef=polyfit(tn,X,1);
LongTermDrift=coef(1)*tn+coef(2);
X=X-LongTermDrift;

if nargin<2
dt=0.001;
end

if nargin>2
	for pf=1:length(varargin)
		if(skip)
			skip=0;
			continue;
		end
		arg=varargin{pf};
		if(ischar(arg))
			arg=lower(arg);
			switch (arg)
			 case {"hann","hn","hp","hm","rect","none","box","ft","bm","exp","kb","pcos","ks","tr","norm","gauss","cgauss","cnorm","nt","pt"}
			  Win_Spec=arg;
			  if(pf<length(varargin))
			  	if(~ischar(varargin{pf+1}))
			  		alpha=varargin{pf+1};
			  		alpha_specified=1;
			  		skip=1;
			  		continue;
			  	end
			  end
			 case {"hanning","hamming","rectangular","flattop","blackman","exponential","poisson","triangular","bartlett","bohman","hann-poisson","kaiser","kaiser-bessel","gaussian","welch","nuttall","planck-taper"}
			  Win_Spec=arg;
			  if(pf<length(varargin))
			  	if(~ischar(varargin{pf+1}))
			  		alpha=varargin{pf+1};
			  		alpha_specified=1;
			  		skip=1;
			  		continue;
			  	end
			  end
			 case {"smooth","sm","pinghua","guanghua"}
			  smooth=1;
			  if(pf<length(varargin))
			  	if(~ischar(varargin{pf+1}))
			  		smpara=varargin{pf+1};
			  		smpara_specified=1;
			  		skip=1;
			  		continue;
			  	end
			  end			  	
			 otherwise
			  plot_option=arg;
			  quickview=1;
			end
		else
			if (~Twin_specified)
				T_winWidth=arg;
				Twin_specified=1;
				continue;
			end
			if (~overlap_specified)
				overlap=arg;
				overlap_specified=1;
				continue;
			end
			
		end
	end
end

if T_winWidth>(L-1)*dt
printf("Window function timescale out of range, default value used (max).\n");
T_winWidth=(L-1)*dt
end

if T_winWidth<0.01*(L-1)*dt
printf("Window function timescale out of range, default value used (min).\n");
T_winWidth=0.01*(L-1)*dt
end

if (overlap>0.9 || overlap<0)
printf("Invalid overlap, please specify from 0.0 to 0.9, default value used.\n");
overlap=0
end

N_win=1+floor(T_winWidth/dt);
n_win=(0:N_win-1);
if size(X,1)>1
n_win=n_win';
end

switch (Win_Spec)
 case {"hanning", "hann", "hn"}
	Win_H=0.5-0.5*cos(2*pi/(N_win-1)*n_win);%Win_recov=2;
 case {"rectangular","rect","none","box"}
	Win_H=ones(size(n_win));%Win_recov=1;
 case {"hamming", "hm"}
	Win_H=0.54-0.46*cos(2*pi/(N_win-1)*n_win);%Win_recov=1.852;
 case {"flattop", "ft"}
	Win_H=0.21557895-0.41663158*cos(2*pi/(N_win-1)*n_win)+0.277263158*cos(4*pi/(N_win-1)*n_win)-0.083578947*cos(6*pi/(N_win-1)*n_win)+0.006947368*cos(8*pi/(N_win-1)*n_win);%Win_recov=4.55;
 case {"blackman", "bm"}
	Win_H=7938/18608-8240/18608*cos(2*pi/(N_win-1)*n_win)+1430/18608*cos(4*pi/(N_win-1)*n_win);%Win_recov=2.381;
 case {"exp", "exponential", "poisson"}
	Win_H=exp(-abs(n_win-(N_win-1)/2)/(N_win-1));%Win_recov=1.2724;
 case {"kaiser", "ks"}
 	if (~alpha_specified)
 	 fprintf("This window function needs to specify alpha.\nDefault value has been used\n");
 	 alpha=10
 	end
	Win_H=besseli(0,alpha*pi*sqrt(1-(2*n_win/(N_win-1)-1).^2))/besseli(0,alpha*pi);
 case {"triangular", "bartlett", "tr"}
	if  ~mod(N_win, 2)
	N_win=N_win-1;n_win(end)=[];
	end
	Win_H=1-2*abs(n_win-(N_win-1)/2)/(N_win-1);%Win_recov=2.0;

 case {"bohman"}
	if  ~mod(N_win, 2)
	N_win=N_win-1;n_win(end)=[];
	end
	Win_H=(1-abs(2*n_win/(N_win-1)-1)).*cos(pi*abs(2*n_win/(N_win-1)-1))+1/pi*sin(pi*abs(2*n_win/(N_win-1)-1));%Win_recov=2.5;

 case {"hann-poisson", "hp"}
	if  ~mod(N_win, 2)
	N_win=N_win-1;n_win(end)=[];
	end
	alpha=0.5;
	Win_H=0.5*(1-cos(2*pi*n_win/(N_win-1))).*exp(-2*alpha/(N_win-1).*abs(n_win-(N_win-1)/2));%Win_recov=2.3256;

 case {"kaiser-bessel", "kb"}
	Win_H=0.402-0.498*cos(2*pi*n_win/(N_win-1))+0.098*cos(4*pi*n_win/(N_win-1))-0.001*cos(6*pi*n_win/(N_win));%Win_recov=2.5;

 case {"pcos"}
 	if (~alpha_specified)
 	 fprintf("This window function needs to specify alpha.\nDefault value has been used\n");
 	 alpha=3
 	end
	Win_H=cos(pi*n_win/(N_win-1)-pi/2).^alpha;%Win_recov=2.381;

 case {"gaussian", "norm", "gauss"}
 	if (~alpha_specified)
 	 fprintf("This window function needs to specify alpha.\nDefault value has been used\n");
 	 alpha=10/3
 	end
	sigma=1/alpha;
	Win_H=exp(-0.5*((2*n_win/(N_win-1)-1)/sigma).^2);%Win_recov=2.703;

 case {"cgauss", "cnorm"}
 	if (~alpha_specified)
 	 fprintf("This window function needs to specify alpha.\nDefault value has been used\n");
 	 alpha=2
 	end
	sigma=1/alpha;
	Win_H=exp(-0.5*((2*n_win/(N_win-1)-1)/sigma).^2)-exp(-0.5*((-1/(N_win-1)-1)/sigma).^2)/(exp(-0.5*((2*(N_win-0.5)/(N_win-1)-1)/sigma).^2)+exp(-0.5*((2*(-N_win-0.5)/(N_win-1)-1)/sigma).^2))*(exp(-0.5*((2*(n_win+N_win)/(N_win-1)-1)/sigma).^2)+exp(-0.5*((2*(n_win-N_win)/(N_win-1)-1)/sigma).^2));%Win_recov=1.7544;

 case {"welch"}
 	Win_H=1-(2*n_win-N_win+1).^2/N_win^2;
 
 case {"nuttall", "nt"}
 	Win_H=0.355768-0.487396*cos(2*pi*n_win/(N_win-1))+0.144232*cos(4*pi*n_win/(N_win-1))-0.012604*cos(6*pi*n_win/(N_win-1));
 
 case {"planck-taper", "pt"}
	if (~alpha_specified)
 	 fprintf("This window function needs to specify alpha.\nDefault value has been used\n");
 	 alpha=2.5
 	end
	sigma=1/alpha;
	ka=(1:floor(sigma*(N_win-1))-1);kb=(ceil((1-sigma)*(N_win-1))+1:N_win-2);
	Win_H=[0,1./(1+exp(sigma*(N_win-1)*(1./ka+1./(ka-sigma*(N_win-1))))),ones(1,length((floor(sigma*(N_win-1)):ceil((1-sigma)*(N_win-1))))),1./(1+exp(sigma*(N_win-1)*(1./(N_win-1-kb)+1./((1-sigma)*(N_win-1)-kb)))),0];

 otherwise
	printf("Unknown Window Type, re-directed to Hanning.\n");
	Win_H=0.5-0.5*cos(2*pi/(N_win-1)*n_win);Win_recov=2;
end

Win_recov=length(Win_H)/sum(Win_H);

CMFFT=zeros(size(n_win));
step_win=floor((1-overlap)*N_win);
count=length([N_win:step_win:L]);
for pf=1:count
CMFFT=CMFFT+fft(Win_H.*X(1+(pf-1)*step_win:N_win+(pf-1)*step_win));
end

ASD=2*Win_recov*abs(CMFFT(1:N_win/2+1))*sqrt(dt/(N_win)/2)/count;
ASD=ASD(2:end-1);

freq=1/dt*(0:(N_win/2))/N_win;
freq=freq(2:end-1);
if size(X,1)>1
freq=freq';
end

if smooth
	fprintf("Smoothing mode activated, please ensure your data longer than size of 100.\n");
	freqs=zeros(size(freq));ASDs=zeros(size(ASD));nReserved=99;
	if smpara_specified
		if length(smpara)>1
			nReserved=smpara(1);
			nLogar=smpara(2);
		else
			nLogar=smpara;
		end
	else
		nLogar=900;
	end
	nidx=[(1:nReserved),floor(logspace(log10(nReserved+1),log10(length(freq)),nLogar+1))];
	nsid=1;
	for pf=2:length(nidx)-1
		if nidx(pf-1)+2>nidx(pf+1)
			continue;
		end
		nSegm=(nidx(pf-1)+1:nidx(pf+1)-1);
		freqs(nsid)=sum(freq(nSegm))/length(nSegm);
		ASDs(nsid)=sum(ASD(nSegm))/length(nSegm);
		nsid=nsid+1;
	end
	freq=freqs(1:nsid-1);
	ASD=ASDs(1:nsid-1);
	clear freqs;
	clear ASDs;
	divfreq=freq(nReserved+1);
	fprintf("Linear Segment: %d\nLogari Segment: %d\ndivided at freq: %e Hz\n\n",nReserved,nsid-nReserved-1,divfreq);
end

if quickview
plot(freq,ASD,plot_option);
set(gca,'Box','on','xminortick','on','yminortick','on','TickDir','out','TickLength',[.02 0]);set(gca,'Xscale','log','Yscale','log','LineWidth',3,'fontsize',16,'fontweight','bold');xlabel('frequency (Hz)');ylabel('Amplitude (X{\cdot}Hz^{-1/2})');
end

end
