%Critical notes:
%.mat file MUST be  named "date_xxx.mat" where xxx are 3 numbers (001,
%012, etc.), y is a single number(1,2,3) and date  is formatted as mmddyy (062421)
%EXAMPLE: "062421_003.mat"


%Load file and read the relevant variable names from .mat file
[fileName, pathName] = uigetfile('*.mat');
fileNameT = fileName(1:end-4);
fileNameSave =  str2num(fileName(8:10));
load(fileName)
a = matfile(fileName);
whoa = who('Ch*');
whob = (whoa{1});
day = str2num(fileName(1:6));

%create dialog for user to input the first template channel to be used from
%Spike2 file
% prompt = {'Enter air pulse length in  seconds', 'Enter drug information (1 for control, 2 for drug)'};
% dlgtitle = 'Spike2 data input';
% startch = inputdlg(prompt,dlgtitle);


totchan = length(who('Ch*')); %detect how many channels are loaded in from Spike2
% airtime = str2num(startch{1})+0.5; %length of air pulse, used to index from gap-free start (0.5 is delay)
% drug = str2num(startch{2});

%Hardcoded airtime and drug for ease of use; be sure to change if either
%condition changes later on
airtime = 5.5;
drug =1;

%Put the Spike2 matlab structure inside a single structure so it can be referenced later on
store = struct;
for x = 1:totchan
    c = a.(['Ch' num2str(x)]);
    store.(['Ch' num2str(x)]) = c;
end
sr = 1/(store.Ch1.interval); %sampling rate in Hz, taken from Spike2 file

%find the first template channel, bath channel, and gap-free for air pulse recordings
timechan = [];
bathchan = [];
for t = 1:totchan
    tempchan(t) = contains(store.(['Ch' num2str(t)]).title, "nw-");
    bathchancheck(t) = double(store.(['Ch' num2str(t)]).title == "Bath Temp");
    airchan(t) = contains(store.(['Ch' num2str(t)]).title, "Gap");
end
a = find(tempchan == 1);
startchan = a(2);
bathchan = find(bathchancheck ==  1, 1, 'first');
airchan  = find(airchan ==1, 1, 'first');
gapcheck = length(airchan);
channelda = store.(['Ch' num2str(startchan)]).title;
channeldata = channelda(4);

holdtable = [];
data = str2double(channeldata);
datachan = (['Ch' num2str(data)]);

if  gapcheck == 0  | length(store.(['Ch' num2str(airchan)]).times) > 50
    spont = 1;
else
    spont =  0;
end
tablesavepath = pwd;
% if spont == 0
%     tablesavepath = "/Volumes/My Passport/Skin_Prep/AIR_EVOKED/";
%     mkdir("/Volumes/My Passport/Skin_Prep/AIR_EVOKED")
% else
%     tablesavepath = "/Volumes/My Passport/Skin_Prep/SPONTANEOUS/";
%     mkdir("/Volumes/My Passport/Skin_Prep/SPONTANEOUS")
% end

if spont == 0
    puffstart = store.(['Ch' num2str(airchan)]).times;
    puffend = store.(['Ch' num2str(airchan)]).times + airtime + 0.5;
else
    puffstart = 1;
    puffend = store.(['Ch' num2str(data)]).length/sr;
end



holdcell = {};
for q  = startchan:totchan
    p = q - (startchan-1);
    chname =  (['Channel' num2str(q)]);
    times = store.(['Ch' num2str(q)]).times;
    values = store.(['Ch' num2str(q)]).values;
    bathval = store.(['Ch' num2str(bathchan)]).values;
    valuesx = length(values(1,:));
    channel = store.([datachan]).values;
    spiketitle = string(p);
    t1 = times; %the vector with the time found by the broken out waveform channel with single spikes
    tt1 = t1*sr; %converting the time vector into indices
%     spikeplot = figure;
    hold on;
    d = values;
    dst = times;
    cols = size(values,2);
    thresh = 10;
    clear peak pk amp
    pk = [];

    %setup for finding spike width
    riseslope = [];
    risetime = [];
    spikemat = [];
    testmat = [];
    slope = [];
    holdw = [];
    area = [];
    searchmat  = [];
    newval = [];
    for i = 1:size(d,1)
        first = channel(round(tt1(i)));
        offset = 0-first;
        newval = channel+offset;
        spikemat(i,:) = newval(round((tt1(i)):round(tt1(i)+valuesx)));%make spikemat the new values
%         spikeplot = plot(newval(round((tt1(i)):round(tt1(i)+valuesx)))); %plotting all the same spikes from  a single broken out channel
        maxi = max(spikemat(i,:));
        mini = min(spikemat(i,:));
        left =  (length(d)/6)*2;
        right = (length(d)/6)*4;
        if abs(maxi) > abs(mini)
            extremetime = find(spikemat(i,:) ==  max(spikemat(i,:)), 1, 'first');
            [maxs maxidx] = max(spikemat(i,:));
        else
            extremetime = find(spikemat(i,:) ==  min(spikemat(i,:)), 1, 'first');
            [maxs maxidx] = min(spikemat(i,:));
        end

        searchmat = [];
        %slope and rise time calculations (10% to 90%)
        bottom = 0.1*(maxs);
        top = 0.9*(maxs);
        endsearch = maxidx;
        searchmat(i,:) = spikemat(i,1:endsearch);
        lefttime = find(abs(searchmat(i,:) - bottom) == min(abs(searchmat(i,:) - bottom)), 1, 'first');
        righttime = find(abs(searchmat(i,:) - top) == min(abs(searchmat(i,:) - top)), 1, 'first');
        xdiff  = righttime - lefttime;
        ydiff = abs(spikemat(i,righttime) - spikemat(i,lefttime));
        riseslope(i) = ydiff/xdiff;
        risetime(i) = ((righttime - lefttime)/sr)*1000;

        %area calculation  (using average of first 5 and last 5 points as a baseline)
        roffset = round(size(d,2))/5;
        if maxidx+round(roffset) < size(spikemat,2)
            firstlast = mean(spikemat(i, [1:5 end-5:end]));
            [lbaseval, lbaseidx] = min(abs(spikemat(i,1:maxidx) - firstlast));
            roffset = round(size(d,2))/5;
            [rbaseval, rbaseidx] = min(abs(spikemat(i, maxidx:maxidx+round(roffset)) - firstlast));
            trueright = rbaseidx + maxidx;
            %plot1 = plot(trueright, spikemat(i,trueright), 'r*')
            %plot2 = plot(lbaseidx, spikemat(i, lbaseidx), 'b*')
            area(i) = sum(abs(spikemat(i,lbaseidx:trueright) - firstlast));
        else
            area(i) = 0;
        end
        rtindex(i) = trueright;
        ltindex(i) = lbaseidx;

        pk(i) = maxidx;
        peak(i) = (pk(i)/sr)+dst(i);
        amp(i) = maxs;
        %     spkwidth =

        if mean(slope) > 0
            [pks,locs,widths,proms] = findpeaks(spikemat(i,:),'MinPeakDistance',valuesx-1);
            if not(isempty(locs))
                ws = widths;
                holdw(i) = ws;
                widthsm = (holdw/sr)*1000;
            end
        else
            [pks,locs,widths,proms] = findpeaks(-1*(spikemat(i,:)),'MinPeakDistance',valuesx-1);
            if not(isempty(locs))
                ws = widths;
                holdw(i) = ws;
                widthsm = (holdw/sr)*1000;
            end
        end
        widths = (holdw/sr)*1000;

        drect = abs(d(i,:)); ... rectifying the voltage trace
            drect = drect(floor(cols/3):floor(2*cols/3)); ... taking the middle 1/3 where the spike should be
            drect = drect(drect>thresh); ... isolating the signal to the portion outside the noise, i.e. above threshold

        spkarea(i) = sum(drect);... a pseudo area calculation "outside" the noise level, not a true area

    end

    pinkpk = peak;
    pamp = amp;
    pwidth = widths;
    pslope = riseslope;
    prise = risetime;
    parea = spkarea;
    pareatrue = area;
%     figure;
%     hold on;
%     pulseplot = plot(pinkpk,[0 1./diff(pinkpk)],'m*');

    dpk = pinkpk;
    dpkdiff = [0 1./diff(pinkpk)];
    damp = pamp;
    darea = parea;
    dareatrue = pareatrue;
    dwidth = pwidth;
    drise = risetime;
    dslope = pslope;
    timesp = t1;
    puffmat = [];
    for i = 1:size(puffstart,1)
        pstart = puffstart(i);
        if spont == 0
            pend = pstart+airtime;
        else
            pend = store.(['Ch' num2str(data)]).length/sr;
        end

        puffd = dpk(dpk>pstart & dpk<pend);

        spikenumhold = min(find((dpk>pstart & dpk<pend)>0));

        %         if i == size(puffstart,1)
        %             range = dpk>pstart & dpk<pend;
        %             holdpkdiff = dpkdiff(range(1:end-1));
        %         else
        %             holdpkdiff = dpkdiff(dpk>pstart & dpk<pend);
        %         end
        holdpkdiff = dpkdiff;
        holddamp = damp(dpk>pstart & dpk<pend);
        holddarea = darea(dpk>pstart & dpk<pend);
        holdwidth = dwidth(dpk>pstart & dpk<pend);
        holdrise = drise(dpk>pstart & dpk<pend);
        holdtime = timesp;
        holdareatrue = dareatrue(dpk>pstart & dpk<pend);
        holdslope = dslope(dpk>pstart & dpk<pend);
        holdarealeft = ltindex(dpk>pstart & dpk<pend);
        holdarearight = rtindex(dpk>pstart & dpk<pend);

        if length(puffd)
            pknum = length(puffd);
            holdpk = puffd;
            for j = 1:pknum
                puffmat(j,1,i) = i; ... Puff Number
                    puffmat(j,2,i) = j+spikenumhold-1; ... absolute Spike Number during whole file
                    puffmat(j,3,i) = holdpk(j); ... spike time in abs seconds (not during the puff)
                    puffmat(j,4,i) = holdpk(j) - pstart; ... spike time from puff initiation
                    %                     if puffd(1) == dpk(1)
                %                     puffmat(j,5,i) = 0; ... inst freq in Hz
                %                     else
                if j == 1
                    puffmat(j,5,i) = 0; ... inst freq in Hz
                else
                    puffmat(j,5,i) = holdpkdiff(j); ... inst freq in Hz
                end
                puffmat(j,6,i) = holddamp(j); ... amplitude in uV
                    puffmat(j,7,i) = holddarea(j); ... spike area
                    puffmat(j,8,i) = holdwidth(j);
                puffmat(j,9,i) = mean(bathval);
                puffmat(j,10,i) = drug;
                puffmat(j,11,i) = holdrise(j);
                puffmat(j,12,i) = spiketitle;
                puffmat(j,13,i) = fileNameSave;
                puffmat(j,14,i) = holdtime(j);
                puffmat(j,15,i) = holdareatrue(j);
                puffmat(j,16,i) = holdslope(j);
                puffmat(j,17,i) = holdarealeft(j);
                puffmat(j,18,i) = holdarearight(j);
                puffmat(j,19,i) = day;
                puffmat(j,20,i) = channeldata;
                puffmat(j,21,i) = p;
                %                 puffmat(j,7) = ; ... spike area


            end
        end
    end

    puffmatwrite = [];
    precell = {};
    for i = 1:size(puffmat,3)
        if i == 1
            holdmat = puffmat(:,:,1);
            puffmatwrite = [(i*ones(size(holdmat,1),1)) holdmat];
        else
            holdmat = [(i*ones(size(holdmat,1),1)) puffmat(:,:,i)];
            puffmatwrite = [puffmatwrite; holdmat];
        end
        precell = puffmatwrite;
    end
    holdcell{1+q-startchan} = precell;
end

titles = {'PuffNum', 'PuffNum2', 'SpikeNum', 'SpikeTime', 'SpikeTimeRel', 'InstFreq', 'SpikeAmp', ...
    'SpikeArea', 'SpikeWidth', 'BathTemp', 'Drug', 'RiseTime', 'SpikeCh', 'File', 'SpikeAbsTime', 'TrueArea', ...
    'RiseSlope', 'AreaIndLeft', 'AreaIndRight', 'Day', 'DataChannel', 'SpikeClass'};
%
% for i = 1:size(holdcell,2)
for j = 1:size(titles,2)
    writedata =  vertcat(holdcell{:});
    %writecell(writedata, ([fileNameT, 'SpikeCh', i]));
    table = array2table([writedata],'VariableNames', {'PuffNum', 'PuffNum2', 'SpikeNum', 'SpikeTime', 'SpikeTimeRel', ...
        'InstFreq', 'SpikeAmp', 'SpikeArea', 'SpikeWidth', 'BathTemp', 'Drug', 'RiseTime', 'SpikeCh', 'File', 'SpikeAbsTime', ...
        'TrueArea', 'RiseSlope','AreaIndLeft', 'AreaIndRight', 'Day', 'DataChannel', 'SpikeClass'});
    tablepath = strjoin([tablesavepath string(fileNameT)]);
    writetable(table, tablepath)
end
% end

display 'Complete'
fileName
clear all
close all