%Try to fit fourier series to motion

%Load data
specStruct.datasetName = 'healthy1';
specStruct.patient = [3];
specStruct.session = [1];
specStruct.exerciseAcceptPrefix = {'KEFO'};
specStruct.exerciseAcceptSuffix = {'SLO'}; % exercise suffix that you want. This will load only exercises that end with 'SLO'
specStruct.exerciseRejectPrefix = {'ARTT','LUNG','GAIT','BKBR','STEP','SDST','STSO','STRA','MEAS','KHEL','KICK','STND','STAI'};  % exercise prefix that you don't want
specStruct.exerciseRejectSuffix = {};  % exercise suffix that you don't wantwant.

options.headerMocap = [];
options.headerImu = [];
%options.headerSegManual = [];
options.headerSegCrop = [];
options.headerSegAlg = [];

if(~exist('patientStruct','var'))
    patientStruct = patientLoader(specStruct, options);
    patientStruct{1}.load();
    pt = patientStruct{1};
    ex = pt.exercises{1};
end
time = ex.dataEkf.time - ex.dataEkf.time(1);

start_times = ex.dataSegManual.timeStart - ex.dataEkf.time(1);
end_times = ex.dataSegManual.timeEnd - ex.dataEkf.time(1);

last_index = find(time < 10);
last_index = last_index(end);
time = time(last_index : end);
Q4 = ex.dataEkf.Q(last_index:end,4);

%Remove frequencies below 0.2 Hz from motion
[a b] = butter(5,[0.1/128*2 10/128*2]);
Q4 = filtfilt(a,b,Q4);

clf
plot(time,Q4);
N = numel([start_times; end_times]);
hold on
plot([start_times; end_times],ones(1,N),'r*');

%Now run through the data and fit it to our 5 frequencies
options = statset('TolFun',5e-3,'MaxIter',30);
window = 25;    %To each side of sample
if 1
    coeffs = zeros(numel(time),12);
    for i =window+1 : numel(time)-window-1

        disp(['Itr ' num2str(i) ' out of ' num2str(numel(time))] );
        dat = Q4(i-window:i+window);
        coeffs(i,:) = nlinfit(linspace(0,(window*2+1)/128,window*2+1),...
            dat,@modelfun,[ones(1,6) zeros(1,6)],options);

    end

    clf;
    plot(time,Q4)
    hold on
    colors = distinguishable_colors(5);
    for i=window+1:1:numel(time)-window-1

        if(mod(1,2))
           col = [1 0 0]; 
        else
           col = [0 1 0]; 
        end
        plot(time(i-window:i+window),modelfun(coeffs(i,:),...
            linspace(0,(window*2+1)/128,window*2+1)),'Color',col);

       for j=1:5
            plot(time(i),mod(coeffs(i,j+5),pi*2),'*','Color',colors(j,:));
       end
    end
end


%Lets try to plot everything in circles

%Seperate points into motion and segments
non_segments = [];
for k=1:numel(start_times)
   indxs = find((time > start_times(k)+0.2) & (time < end_times(k)-0.2));
   non_segments = [non_segments; indxs];
end

for j=1:6
    subplot(3,2,j);
    hold on
    plot3(cos(coeffs(:,j+6)),sin(coeffs(:,j+6)),abs(coeffs(:,j)),'r*');
    plot3(cos(coeffs(non_segments,j+6)),sin(coeffs(non_segments,j+6)),abs(coeffs(non_segments,j)),'b*');
    grid on
    xlabel('x');
    ylabel('y');
    zlabel('z');
    view([20 20]);
end

clf
for i = window+1:3:numel(time)-window-1
    col = [1 0 0];
    for k=1:numel(start_times)
       if(time(i) > start_times(k)+0.2&& time(i) < end_times(k)-0.2);
           col = [0 0 1];
           break;
       end
    end
    
end

