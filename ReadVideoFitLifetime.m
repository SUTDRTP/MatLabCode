% Read videos and measure RTP lifetime
%In this version of the program, we modified the code to fit the decay to
% a*exp(-bt) + c

close all;
clear all;
clc;


maxFrame=5000; %The max numbmer of images to process


%plot setting
lineWidth=3;
markerSize=8;
fontSize=16;

%Specify the video file name
% videoName='SamsungS10liteWR1error.mp4'; %video name

[videoName,path,indx] = uigetfile( ...
    {   '*.*',  'All Files (*.*)'}, ...
    'Select a Video File');


%Other Settings
COI=2; %Channel of Interest: 1 for Red; 2 for Green; 3 for Blue.
selectRegion= 0; %if 0, the entire image is selected; else user can define the region of interest.

%Define Output File Names
dotIdx=find(videoName=='.');
temp=videoName(1:dotIdx-1); %a string to name output file

fileName1=[temp 'RawSignal'];
fileName2=[temp 'RedChannelDecay'];
fileName3=[temp 'GreenChannelDecay'];
fileName4=[temp 'BlueChannelDecay'];

timeOffSet=150; %ignore images taken within 150 ms after switching off the UV lamp


%% Load Video Data

v = VideoReader(videoName);

%Determine the number of images in the video and the framerate
nFrame=v.NumFrames;

%Determine the offset number (number of images to discard after the UV lamp is off
%This offset accounts for sharking or other noises caused while switching
%off the lamp
frameRate=v.FrameRate; %based on camera setting
frameTime=1000/frameRate; %time for taking one image
indexOFFSET=ceil(timeOffSet/frameTime); %images to ignore after switching off UV lamp (considering time delay!)

if selectRegion
    % select the region of interest based on the first; draw a big box
    % considering motions/vibrations
    for i=1:1 %nFrame
        FirstFrame = read(v,i);

        figure;
        imshow(FirstFrame);
        % Rectangle position is given as [xmin, ymin, width, height]
        objectRegion=round(getPosition(imrect))
    end


    % save emission intensity info
    figure;
    for i=1:min(nFrame, maxFrame)

        msg=['Frame #' num2str(i) ' of ' num2str(min(nFrame, maxFrame)) ' images']

        frame = read(v,i);

        % Select part of the image
        img_cropped = frame(objectRegion(2) + (0:objectRegion(4)), objectRegion(1) + (0:objectRegion(3)), :);

        
        imshow(img_cropped);
        

        %save total intensity;
        intensity(1,i)=double(sum(sum(img_cropped(:,:, 1)))); %total intensity in blue channel
        intensity(2,i)=double(sum(sum(img_cropped(:,:, 2)))); %total intensity in yellow channel
        intensity(3,i)=double(sum(sum(img_cropped(:,:, 3)))); %total intensity in red channel

    end

else
    % save emission intensity info
    figure;
    for i=1:min(nFrame, maxFrame)

        msg=['Frame #' num2str(i) ' of ' num2str(min(nFrame, maxFrame)) ' images']

        frame = read(v,i);

        imshow(frame);

        %save total intensity;
        intensity(1,i)=double(sum(sum(frame(:,:, 1)))); %total intensity in blue channel
        intensity(2,i)=double(sum(sum(frame(:,:, 2)))); %total intensity in yellow channel
        intensity(3,i)=double(sum(sum(frame(:,:, 3)))); %total intensity in red channel

    end

end


close all;

% save time info based on frame rate and frame numbers
timeIdx=1:1:min(nFrame, maxFrame);
timeIdx=double(timeIdx);
time=(1000/frameRate)*timeIdx;

'Video data retrieval is completed'
clear FirstFrame
clear frame


%% determine starting time t0 and remove data before that
look4Start=diff(intensity(COI,:)); %perform differentiation
index0=find(look4Start==min(look4Start)); %choose the most negative slop; this intensity drop corresponds to switching off the UV lamp.

index=index0+indexOFFSET; %the indexOFFSET accounts for noises while switching off the lamp (i.e., shaking, time-lag etc.)

useTime=time(index:end);
useTime=useTime-useTime(1);
useIntensity=intensity(:, index:end);

% plot raw intensity data
figure;
hold on;
plot(time, intensity(1,:), 'r-o', 'markersize', markerSize, 'markeredgecolor', 'r', 'markerfacecolor', 'r', 'linewidth', lineWidth);
plot(time, intensity(2,:), 'g-o', 'markersize', markerSize, 'markeredgecolor', 'g', 'markerfacecolor', 'g', 'linewidth', lineWidth);
plot(time, intensity(3,:), 'b-o', 'markersize', markerSize, 'markeredgecolor', 'b', 'markerfacecolor', 'b', 'linewidth', lineWidth);
plot([index index]*1000/frameRate, [0 max(max(intensity))*1.1], 'k--', 'linewidth', lineWidth);
legend('Red Channel', 'Green Channel', 'Blue Channel', 'Starting Point');
set(gca, 'fontname', 'Arial', 'fontsize', fontSize);
xlabel('Time (ms)', 'fontname', 'Arial', 'fontsize', fontSize);
ylabel('Emission Intensity', 'fontname', 'Arial', 'fontsize', fontSize);
myTitle='Raw Emission decay dynamics ';
title(myTitle, 'fontname', 'Arial', 'fontsize', fontSize);

% saveas(gcf,[fileName1 '.fig']);
saveas(gcf,[fileName1 '.png']);
saveas(gcf,fileName1, 'epsc');

% plot selected intensity data
figure;
hold on;
plot(useTime, useIntensity(1,:), 'r-o', 'markersize', markerSize, 'markeredgecolor', 'r', 'markerfacecolor', 'r', 'linewidth', lineWidth);
plot(useTime, useIntensity(2,:), 'g-o', 'markersize', markerSize, 'markeredgecolor', 'g', 'markerfacecolor', 'g', 'linewidth', lineWidth);
plot(useTime, useIntensity(3,:), 'b-o', 'markersize', markerSize, 'markeredgecolor', 'b', 'markerfacecolor', 'b', 'linewidth', lineWidth);
legend('Red Channel', 'Green Channel', 'Blue Channel');
set(gca, 'fontname', 'Arial', 'fontsize', fontSize);
xlabel('Time (ms)', 'fontname', 'Arial', 'fontsize', fontSize);
ylabel('Emission Intensity', 'fontname', 'Arial', 'fontsize', fontSize);
myTitle='Emission decay dynamics ';
title(myTitle, 'fontname', 'Arial', 'fontsize', fontSize);


%% perform the red channel fitting
myIntensity=useIntensity(1,:);

% Perform fitting without removing the offset
f = @(b,x) b(1).*exp(b(2).*x) + b(3);
nrmrsd = @(b) norm(myIntensity - f(b,useTime));                          % Residual Norm Cost Function
B0 = rand(3,1);
B0(1)=myIntensity(1);
B0(3)=myIntensity(end);
B0(2)=(0-log(myIntensity(1)-myIntensity(end)))/useTime(end);

% B0(2)=-1/100;
% Choose Appropriate Initial Estimates
[B,rnrm] = fminsearch(nrmrsd, B0);
myIntensityFit=B(1).*exp(B(2).*useTime) + B(3);

tao_red=-1/B(2);
['tao_red = ' num2str(tao_red) ' ms']
myExpIntensity=useIntensity(1,:);
myFitIntensity=myIntensityFit;

figure;
hold on;
plot(useTime, myExpIntensity, 'ro', 'markersize', markerSize, 'linewidth', lineWidth);
plot(useTime, myFitIntensity, 'k-', 'linewidth', lineWidth );
legend('experimental', 'fitted');
set(gca, 'fontname', 'Arial', 'fontsize', fontSize);
xlabel('Time (ms)', 'fontname', 'Arial', 'fontsize', fontSize);
ylabel('Emission Intensity', 'fontname', 'Arial', 'fontsize', fontSize);
myTitle='Decay dynamics in the red channel';
title(myTitle, 'fontname', 'Arial', 'fontsize', fontSize);

% saveas(gcf,[fileName2 '.fig']);
saveas(gcf,[fileName2 '.png']);
saveas(gcf,fileName2, 'epsc');

%% perform the green channel fitting
myIntensity=useIntensity(2,:);

% Perform fitting without removing the offset
f = @(b,x) b(1).*exp(b(2).*x) + b(3);
nrmrsd = @(b) norm(myIntensity - f(b,useTime));                          % Residual Norm Cost Function
B0 = rand(3,1);
B0(1)=myIntensity(1);
B0(3)=myIntensity(end);
B0(2)=(0-log(myIntensity(1)-myIntensity(end)))/useTime(end);

% Choose Appropriate Initial Estimates
[B,rnrm] = fminsearch(nrmrsd, B0);
myIntensityFit=B(1).*exp(B(2).*useTime) + B(3);

tao_green=-1/B(2);
['tao_green = ' num2str(tao_green) ' ms']
myExpIntensity=useIntensity(2,:);
myFitIntensity=myIntensityFit;

figure;
hold on;
plot(useTime, myExpIntensity, 'go', 'markersize', markerSize, 'linewidth', lineWidth);
plot(useTime, myFitIntensity, 'k-', 'linewidth', lineWidth );
legend('experimental', 'fitted');
set(gca, 'fontname', 'Arial', 'fontsize', fontSize);
xlabel('Time (ms)', 'fontname', 'Arial', 'fontsize', fontSize);
ylabel('Emission Intensity', 'fontname', 'Arial', 'fontsize', fontSize);
myTitle='Decay dynamics in the green channel';
title(myTitle, 'fontname', 'Arial', 'fontsize', fontSize);

% saveas(gcf,[fileName3 '.fig']);
saveas(gcf,[fileName3 '.png']);
saveas(gcf,fileName3, 'epsc');

%% perform the blue channel fitting
myIntensity=useIntensity(3,:);

% Perform fitting without removing the offset
f = @(b,x) b(1).*exp(b(2).*x) + b(3);
nrmrsd = @(b) norm(myIntensity - f(b,useTime));                          % Residual Norm Cost Function
B0 = rand(3,1);
B0(1)=myIntensity(1);
B0(3)=myIntensity(end);
B0(2)=(0-log(myIntensity(1)-myIntensity(end)))/useTime(end);

% Choose Appropriate Initial Estimates
[B,rnrm] = fminsearch(nrmrsd, B0);
myIntensityFit=B(1).*exp(B(2).*useTime) + B(3);

tao_blue=-1/B(2);
['tao_blue = ' num2str(tao_blue) ' ms']
myExpIntensity=useIntensity(3,:);
myFitIntensity=myIntensityFit;

figure;
hold on;
plot(useTime, myExpIntensity, 'bo', 'markersize', markerSize, 'linewidth', lineWidth);
plot(useTime, myFitIntensity, 'k-', 'linewidth', lineWidth );
legend('experimental', 'fitted');
set(gca, 'fontname', 'Arial', 'fontsize', fontSize);
xlabel('Time (ms)', 'fontname', 'Arial', 'fontsize', fontSize);
ylabel('Emission Intensity', 'fontname', 'Arial', 'fontsize', fontSize);
myTitle='Decay dynamics in the blue channel';
title(myTitle, 'fontname', 'Arial', 'fontsize', fontSize);

% saveas(gcf,[fileName4 '.fig']);
saveas(gcf,[fileName4 '.png']);
saveas(gcf,fileName4, 'epsc');

%% Save Data for future analysis
save(temp);
