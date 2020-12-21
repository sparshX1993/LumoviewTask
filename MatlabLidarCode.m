%% LIDAR 2D DATA REPRESENTATION : Author: Sparsh Tiwari : Date : 21/12/2020
% Reading the data from the text file and distributing it in various
% channels and getting out good quality data based on light quality
%%
clear all,close all
fileID = fopen('out_startplatz_cut.txt');
C = textscan(fileID,'%s %d %f %f');
NewRotation=C{1,1};
MeasurementQualities=C{1,2};
MeasurementAngles=C{1,3};
MeasurementDistance=C{1,4};
Size=size(MeasurementDistance);%getting the size of data channels
MeasurementDistanceX=zeros(Size(1),1);%Creating a zero matrix with size of input data
MeasurementDistanceY=zeros(Size(1),1);
k=zeros(Size(1),1);
k1=zeros(Size(1),1);

P=-2*pi; %taking P as constant value of 2*Pi
for i=1:Size(1)
    if (MeasurementQualities(i,1) ==15 ) & (MeasurementDistance(i,1) < 5000) & (MeasurementDistance(i,1) > 1500)
    [MeasurementDistanceX(i,1), MeasurementDistanceY(i,1)] = pol2cart(MeasurementDistance(i,1),P*(MeasurementAngles(i,1))/360);
    else 
        MeasurementDistanceX(i,1)=NaN;
        MeasurementDistanceY(i,1)=NaN;    
        MeasurementDistance(i,1)=NaN;
    end

end
%% Using fillmissing function om the NaN entries of MeasurementDistanceX and MeasurementDistanceY

 X2 =  fillmissing(MeasurementDistanceX,'linear');
 Y2 =  fillmissing(MeasurementDistanceY,'linear');
 
x=MeasurementDistanceX;% moving the value of measurementdistanceX to new variable x
y=MeasurementDistanceY;

x1=find(isnan(x)); %finding out the NaN entries on x and y
y1=find(isnan(y));

x2=zeros(Size(1),1);%Creating a zero matrix with size of input data
y2=zeros(Size(1),1);

x2(x1,1)=X2(x1,1);%Moving the missing data from the linear model values to the new empty matrix
y2(y1,1)=Y2(y1,1);% Using Indexing only NaN values will be replaced with linear values. 

x2(x2==0)=NaN; %changing every zero values to NaN
y2(y2==0)=NaN;
%% Deleting every NaN values from x and y(parameters inheriting the values of measurementDistanceX/Y

x(find(isnan(x)))=[];
y(find(isnan(y)))=[];

%% Plotting for MeasurementDistance X/Y(with missing values) and the x2/y2

plot(MeasurementDistanceX,MeasurementDistanceY,'b-',x2,y2,'r-')
%% Plotting for one single rotation(0-360).In this case ploting for frames 167-334
%% we can also do an average for all the MeasurementDistanceX/Y based on the angle rotation from 0-360 
%% Condensing  the data set for more accurate resutlt. 
plot(x(167:334,1),y(167:334,1),'b-',x2(167:334,1),y2(167:334,1),'r-')
xlabel('Measurement Distance X')
ylabel('Measurement Distance Y')
title('2D Lidar Data representation')
legend('Original Data','Reconstructed Data','Location','southeast')
%% Exporting Original Measurement Data with New Measurement Data into a CSV file
filename1=[MeasurementDistanceX X2 MeasurementDistanceY Y2]
writematrix(filename1, 'MeasurementData.csv');
