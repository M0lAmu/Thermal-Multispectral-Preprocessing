%%%% INPUT %%%
LocIm = 'G:\Molamu processing\202324 Season\Thermal\NOvember\22 Nov\tif\'; %location of the images

%v2: all images are exported from Thermoviewer, the sharpest per second
%is/are selected; the other images are removed (to save space)
%(Opt_RemUnsh)

%v3: v2 + also working when no images with GPS coordinates are present in the dataset. 

% Options:

Opt_Tcorr = 1;  %put 1 if you want to perform this option. If weather conditions during the flight were highly changing, 
                % use this to correct for Temperature fluctuations which
                % are stored in an additional folder 
Opt_RemUnsh = 1; %Unsharp images are removed in the first step of the process 
Opt_Sharp = 1;  %If you want to move unsharp images to separate folder
Opt_PosMulti = 0; %if you want to estimate the XYZ and PitchYawRoll for the thermal cameras based on the multispectral data. 
                % This requires that the multispectral data are processed
                % in METASHAPE, and that you have exported the estimated
                % camera locations. Matching will be done based on the GNSS
                % time. This works only if you have selected the option to
                % create a XXX_meta.csv file when exporting the cameras in
                % Thermoviewer! Regardless, another file is also created
                % with the camera positions recorded by the GPS.
Opt_LocOnly = 0; %1 if you want to use the script only to save the location of the cameras; otherwise, 1
Opt_Csv = 0;    % 1 for csv-data, 0 for Tif data
               
                
% Parameters:

RH = 57.32;         %Relative Humidity, in percentage. If unknown, try looking it up at https://wow.meteo.be/nl/. Else, put 0
Tair = 22;       %Air temperature, in  C; If unknown, try looking it up at https://wow.meteo.be/nl/. Else, put 0
Tbg = -11;        %Background temperature. Seen in Thermoviewer as Tbr of alu panels. 
FH = 33;        %Flight Height - actually height to object (flight height minus crop/tree height)
ImPerSec = 2;   %Images per second you want to keep.

% Only needed for Opt_PosMulti, ie if you want to estimate the XYZ and
% YawPitchRoll of the thermal based on the multispectral data. This works
% only for the TIF option!
%LocMulti = 'D:\Forbio\20210721_Zedelgem\multi\Samen';       % Location of multispectral Tif data. Careful, this needs to be the original data! (Not the reflectance product)
%LocTxtMulti = 'D:\Forbio\20210721_Zedelgem\Agisoft'; %Location of file with multispectral locationHere
%NameTxt = '20210721_Zed_Multi_Pos.txt';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extra INPUT: don't change unless absolutely necessary and you know what you're doing!
TLimLow = 253; TLimHigh = 353; %in K; Don't change unless absolutely necessary
emiss= 0.985;    %Emissivity

SmoothLength = 60; %Total smoothing for Opt_Tcorr; this is in number of images, with usually 3 images per second. (60= smoothing of 10s before and 10s after point) 

PTileThresh = 3; % Parameters for removing the unsharp data
MaxGap = 4;
StepThresh = 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%This script reads the data from Thermoviewer (either csv or tif, to be mentioned), performs the full radiometric correction, 
%and saves the Ts-data as Tif files. 

%It then looks at sharp vs unsharp images, and moves less sharp images to a
%separate folder. 

%In case Ts-data were obtained in cloudy/changing conditions, a correction
%(Opt_Tcorr) is added, which performs a smoothing of the Ts-data and
%corrects for temperature fluctuations. This has the danger of also
%eliminating areas with hotter/cooler conditions. A smoothing of 60 frames
%(roughly 20 seconds of flight, or about 60-80m) is now implemented. This
%can be adjusted 

%%%%%%%%%%%%%%% Actual Script
warning ('off','all');

%Create LocExp
if ~(strcmp(LocIm(end),'\') | strcmp(LocIm(end),'/'))
    LocIm = [LocIm '\'];
end
%if ~(strcmp(LocTxtMulti(end),'\') | strcmp(LocTxtMulti(end),'/'))
    %LocTxtMulti = [LocTxtMulti '\'];
%end
%if ~(strcmp(LocMulti(end),'\') | strcmp(LocMulti(end),'/'))
    %LocMulti = [LocMulti '\'];
%end
LocExp = [LocIm 'Tif_Ts\']; mkdir(LocExp);

if Opt_Csv == 1
    DirIm = dir([LocIm '*.csv']);
else
    DirIm = dir([LocIm '*.tif']);
end

%% Import some basic data (time, location) per image; define if NoGPS
DirMeta = dir([LocIm '*_meta.csv']);
%if round(length(DirIm)/10,0) == round(length(DirMeta)/10,0) %If by accident you saved a metadata file for each image, this part corrects for this.
   %for k1=1:length(DirMeta)
       %if k1==1 
           %ThermalInf = ImportMeta([LocIm DirMeta(1,1).name]);
       %else
           %ThermalInf = [ThermalInf; ImportMeta([LocIm DirMeta(k1,1).name])];
       %end
   %end
%else
    ThermalInf = ImportMeta([LocIm DirMeta(1,1).name]);
%end
  
 % Check: if no GPS recordings
if ~(any(isnumeric(ThermalInf.Latitude)))
    disp('No GPS data recorded');
    NoGPS=1;
else  
    NoGPS=0;
end


%% First loop: only keep the sharpest X images per second (ImPerSec), remove all other
%images

TTime = nan(length(DirIm), 1);
TSharp = TTime;
TTime_corr = nan(size(TTime));
KeepT = zeros(size(TTime));

if Opt_LocOnly == 0
    for i1 = 1:length(DirIm)
        A = f_readIm([LocIm DirIm(i1, 1).name], Opt_Csv);
        A = double(A);
        [TSharp(i1)]=estimate_sharpness(A);
    end
end
%clear A


%% First, get the time of the thermal data capture
if length(DirIm) ~= length(ThermalInf.Date) 
    
    disp('Nr of Tifs not the same as in Meta-table');
end
%% 
    if NoGPS == 0
                    % write down thermal time
                     %Thermal Time
       DateS = cellstr(ThermalInf.Date(2));
       Date = datetime(str2num(DateS{1}(7:10)), str2num(DateS{1}(4:5)), str2num(DateS{1}(1:2)));
       Time = datetime(ThermalInf.Time, 'InputFormat', 'HH:mm:ss');


                     % Calculate the combined datetime
       TTime = Date + timeofday(Time);
    end

 % Find unique timestamps
 if NoGPS == 0
    [~, Id] = unique(TTime);
 else
     Id = 1:length(DirIm);
 end

 % Keep the sharpest images per second
 for i = 1:length(Id)
     % find indices corresponding to the current second
     indices_in_second = find(TTime >= TTime(Id(i)) & TTime < TTime(Id(i)) + seconds(1));

     % determine the number of images to keep for this second
     num_images_to_keep = min(length(indices_in_second), ImPerSec);

     % find the indices of the sharpest images in this second
     [~, KeepId] = maxk(TSharp(indices_in_second), num_images_to_keep);

     % Mark the selected images for keeping
     KeepT(indices_in_second(KeepId)) = 1;
 end

 % Remove unselected images
 Rem = DirIm(KeepT == 0);
 DirIm = DirIm(KeepT == 1);
                     
            

         % Now remove the unsharp images and continue working with the sharp images only.

         TTime_corr = TTime_corr(KeepT == 1);
         TTime = TTime(KeepT == 1);
         TSharp=TSharp(KeepT == 1);
         ThermalInf = ThermalInf(KeepT == 1,:);
    
         %if Opt_RemUnsh
             %for i2=1:length(Rem)
                 %delete([LocIm Rem(i2).name]);
             %end
         %end
     
    %else
         %If no GPS (or time) stamp, we will group the images per five, and keep the ImPerSec nr of images per group of 5.

           %if Opt_RemUnsh
               %LDir = round(length(DirIm)/5);
               %for i2=1:LDir
                    %IDs = [(i2-1)*5+1:1:i2*5];
                    %[~,KeepId] = maxk(TSharp(IDs),ImPerSec);
                    %for i3=1:ImPerSec
                        %KeepT((i2-1)*5+KeepId(i3))=1;
                    %end
               %end
%end
      
        % Now remove the unsharp images and continue working with the sharp images only.
       
         %TSharp=TSharp(KeepT==1);
         %ThermalInf = ThermalInf(KeepT==1,:);
         %Rem = DirIm(KeepT==0);
         %DirIm = DirIm(KeepT==1);
         
         %if Opt_RemUnsh
             %for i2=1:length(Rem)
             %    delete([LocIm Rem(i2).name]);
             %end
         %end


    %end

%end
%% 
 
 %Convert Tbg and ImCsv to K if needed
  if Tbg <200
     Tbg = 273.15+Tbg;
  end

  if Opt_Tcorr == 1
   TAv = nan(size(DirIm));
  end

% calculating transmission
% Code % Obtained from Heinemann et al 2020 RS https://doi.org/10.3390/rs12071075
%Can be placed in the loop, if the FH is made dynamic. 

     if Tair == 0
         A = f_readIm([LocIm DirIm(1,1).name], Opt_Csv);
         Tair = mean(A(:)); 
         if Tair > 200 %convert to K if needed (but shouldn't be!)
             Tair = Tair-273.15;
         end
     end
     
     if RH == 0
         RH = 50;
     end

omega = RH/100*exp(6.8455*10^-7*Tair.^3 -2.7816*10^-4*Tair.^2+6.939*10^-2*Tair+1.5587);
Tau_atm = 1.9*exp(-sqrt(FH)*(0.0066 -0.0023*sqrt(omega))) + (1-1.9)*exp(-sqrt(FH)*(0.0126 -0.0067*sqrt(omega)));

 
 
%%
%Loop through image directory: reading thermal images, converting to Ts and
%exporting

if ~Opt_LocOnly == 1

    for i1=1:length(DirIm)
    %  try
        A = f_readIm([LocIm DirIm(i1,1).name], Opt_Csv);
        %% Convert temperature data to K
         if A(1,1)<200
            A = 273.15+A;
         end

         %% Correcting for atmospheric conditions
           Ts = ((A.^4 - (1-emiss)*Tbg.^4 - (1-Tau_atm).*(Tair+273.15).^4)/(emiss.*Tau_atm)).^0.25;

                   %Convert images to 16-bit tif
           TsConv = (Ts -TLimLow)./(TLimHigh-TLimLow).*65535;
           Image16 = uint16(round(TsConv));
           NameIm = [LocExp DirIm(i1,1).name(1:end-4) '.tif'];
           imwrite(Image16,NameIm,'tif','Compression','none'); 
           clear NameIm Image16 

       if Opt_Tcorr
            TAv(i1)=trimmean(Ts(:),20);
       end


    %  catch 
    %         display(['Issue with image ' DirIm(i1,1).name(1:end-4)]);
    %  end   
    end


    %Temperature correction, if needed
    if Opt_Tcorr

        %% SmoothingFunction 
        % From Trimmean, calculate mean temperature 
        T_Smooth = nan(size(TAv));
        for i1=1:length(T_Smooth)
            T_Smooth(i1) = mean(TAv(max(1,i1-SmoothLength*0.5):min(length(T_Smooth),i1+SmoothLength*0.5)), 'omitnan');
        end
        %Correction
        TCorr = mean(T_Smooth, 'omitnan') - T_Smooth;

        %% Create directory
        LocCorr = [LocExp 'Ts_Corr\']; mkdir(LocCorr);

        %% loop: load photos, and perform the correction
        for i1=1:length(DirIm)
            %Load images
            A = imread([LocExp DirIm(i1,1).name]);A=double(A);
            %Convert back to Ts
            Ts = (TLimHigh-TLimLow).*A./65535 + TLimLow;
            %Convert images to 16-bit tif
            TsConv = (Ts + TCorr(i1) -TLimLow)./(TLimHigh-TLimLow).*65535;
            Image16 = uint16(round(TsConv));
               NameIm = [LocCorr DirIm(i1,1).name];
               imwrite(Image16,NameIm,'tif','Compression','none');
        end
    end

    %Output metadata
    OutText = {'Tlow_Tif'; 'Thigh_Tif'; 'emissivity'; 'TBg'; 'Tair'; 'RelHum';'Tau_atm'; 'FlightHeight'; 'Opt_Sharp'; 'Opt_Tcorr'};
    OutDat = {TLimLow;TLimHigh;emiss;Tbg;Tair;RH;Tau_atm;FH;Opt_Sharp; Opt_Tcorr};
    OutT = table([OutText OutDat]);
    writetable(OutT,[LocExp 'Parameters.xlsx']);
    if Opt_Tcorr
        OutText = {'Tlow_Tif'; 'Thigh_Tif'; 'emissivity'; 'TBg'; 'Tair'; 'RelHum';'Tau_atm'; 'FlightHeight'; 'Opt_Sharp'; 'Opt_Tcorr'; 'SmoothLength'};
        OutDat = {TLimLow;TLimHigh;emiss;Tbg;Tair;RH;Tau_atm;FH;Opt_Sharp; Opt_Tcorr; SmoothLength};
        writetable(OutT,[LocCorr 'Parameters.xlsx']);
    end
    %% Sharpening data

    if Opt_Sharp
        LocUnsharp = [LocExp 'Unsharp\'];mkdir(LocUnsharp);
        if Opt_Tcorr
            LocUnsharpCorr = [LocCorr 'Unsharp\'];mkdir(LocUnsharpCorr);
        end

        % Step 1: Define sharpness of all images

     % DirIm = dir([LocExp '*.tif']);
      if Opt_Csv == 1
          NameTif = struct2cell(DirIm)';NameTif = NameTif(:,1);
          for i3 =1:length(NameTif)
              NameTif{i3,1}=[NameTif{i3,1}(1:end-4) '.tif'];
          end
      else
        NameTif = struct2cell(DirIm)';NameTif = NameTif(:,1);
      end
% 
%         TSharp = nan(length(DirIm), 1);Gaps = TSharp;
% 
%         for i1=1:length(TSharp)
%             %Open image
%             A =imread([LocExp DirIm(i1,1).name(1:end-4) '.tif']);
%             A = double(A);
%             [TSharp(i1)]=estimate_sharpness(A);
%         end
%         clear i1 A

        %% Step 2 (sharpening): Define optimal sharpness

        PTileThresh_Sel=PTileThresh;
        x=0;
        while x==0
                Thresh = prctile(TSharp,PTileThresh_Sel);
                [~,IDEdge] = min(abs(TSharp-Thresh));
                A =imread([LocIm DirIm(IDEdge,1).name(1:end-4) '.tif']);
                A=double(A);LL = min(A(:)); ML =max(A(:));
                Handles.H=figure('Name','Sharpness');  imshow(A,[LL ML]);
                prompt = 'Sharp enough? (1/0)';
                x = input(prompt);
                close(Handles.H);
            if x==0
                PTileThresh_Sel = PTileThresh_Sel + StepThresh;
            end
        end

        ImageDropped = sum(TSharp < Thresh);
        disp(['Selected threshold at percentile ' num2str(PTileThresh_Sel) ' (at sharpness value of ' num2str(Thresh) '); ' num2str(ImageDropped) ' images of total of ' num2str(length(TSharp)) ' dropped']);

        %% Step 3 (sharpening): select images to move
        %Calculate GAPS
        Gaps = zeros(size(TSharp));
        for i1 = ceil(MaxGap*0.5):(length(DirIm)-ceil(MaxGap*0.5))
            if TSharp(i1) < Thresh
                R = 1;
                %Debugging output
                disp(['i1 = ' num2str(i1)]);
                while (i1+R-1) <= length(TSharp) && TSharp(i1+R-1) < Thresh
                    disp(['R = ' num2str(R)]);
                    disp(['Length of TSharp = ' num2str(length(TSharp))]);
                    
                    R = R + 1;
    %%%%%%%%%%%                % Additional check
                    if (i1 + R - 1) > length(TSharp)
                        disp('Breaking out of loop to avoid out-of-bounds access');
                        break;
                    end
                end
                Gaps(i1)=R;
                clear R
            end
        end

        MoveFile = zeros(size(TSharp));
        MoveFile(TSharp < Thresh)=1;
        for i1=ceil(MaxGap*0.5):length(DirIm)-ceil(MaxGap*0.5)
            if Gaps(i1) > MaxGap & (TSharp(i1)==max(TSharp(max(1,i1-floor(MaxGap*0.5)):min(length(DirIm),i1+floor(MaxGap*0.5)))))
                MoveFile(i1)=0;
            end
        end


        %% Step 4 (sharpening): move tifs
        MoveTif = NameTif(MoveFile==1);
        if ~isempty(MoveTif)
            for i2=1:length(MoveTif)
                a =movefile([LocExp MoveTif{i2}],[LocUnsharp MoveTif{i2}]);
                if Opt_Tcorr
                    a =movefile([LocCorr MoveTif{i2}],[LocUnsharpCorr MoveTif{i2}]);
                end
            end
        end
    end
end

%% Part on estimating camera positions based on multispectral data

%if ~NoGPS

    %TTime = nan(length(DirIm),1);

    % Look for the meta-file 

    % DirMeta = dir([LocIm '*_meta.csv']);
    % if round(length(DirIm)/10,0) == round(length(DirMeta)/10,0) %If by accident you saved a metadata file for each image, this part corrects for this.
    %    for k1=1:length(DirMeta)
    %        if k1==1 
    %            ThermalInf = ImportMeta([LocIm DirMeta(1,1).name]);
    %        else
    %            ThermalInf = [ThermalInf; ImportMeta([LocIm DirMeta(k1,1).name])];
    %        end
    %    end
    % else
    %     ThermalInf = ImportMeta([LocIm DirMeta(1,1).name]);
    % end
    %if Opt_PosMulti

        %TimeDiff = 0;
        %Bias=[0,0, 0,0,0,0]; 
        %try
            %% First, get the time of the thermal data capture
            %if length(DirIm) == length(ThermalInf.Date)
                %Thermal Time
    %             DateS=cellstr(ThermalInf.Date(2)); Date = datenum(str2num(DateS{1}(7:10)),str2num(DateS{1}(4:5)),str2num(DateS{1}(1:2)));
    %             TTime = Date + datenum(ThermalInf.Time) - floor(datenum(ThermalInf.Time));
    %                         
    %             TTime_corr = TTime;
    %             %Scan the nr of images taken per sec.
    %             UnTT = unique(TTime);
    %                 for i7=1:length(UnTT)
    %                     Id = find(UnTT(i7)==TTime);
    %                     for i8 = 1:length(Id)
    %                         TTime_corr(Id(i8)) = TTime(Id(i8)) + (i8-1)/(length(Id).*(24.*3600));
    %                     end
    %                 end


                %% Then, get the time of the multispectral data
                 %ListMulti = dir([LocMulti, '*_6.tif']);
                 %if isempty(ListMulti)
                     %ListMulti = dir([LocMulti, '*_1.tif']);
                 %end
                  %ListMulti=struct2cell(ListMulti)';ListMulti=ListMulti(:,1);
                  %TimeMulti = nan(size(ListMulti));
                  %LabelMulti = cell(size(ListMulti));
                  %for i2=1:length(ListMulti)
                                %ImMulti = imfinfo([LocMulti ListMulti{i2}]); DT = ImMulti.DateTime;
                                %TimeMulti(i2) = datenum(str2num(DT(1:4)),str2num(DT(6:7)),str2num(DT(9:10)),...
                                                        %str2num(DT(12:13)),str2num(DT(15:16)),str2num(DT(18:19)));
                               %LabelMulti{i2} = str2num(ListMulti{i2}(5:8));      
                  %end

                  %% Load the camera positions

                  %NameMulti_Pos =Imp_MultiNames([LocTxtMulti NameTxt]);
                  %NameMulti_Pos = cellstr(NameMulti_Pos);
                  %MultiPos = Import_Pos([LocTxtMulti NameTxt]);
                 %% MultiPos.Names = NameMulti_Pos.MultiPos; %clear NameMulti_Pos
                  %MultiPos.X_est= MultiPos.X_est.*10^-6;
                  %MultiPos.Y_est= MultiPos.Y_est.*10^-6;
                  %MultiPos.Z_est= MultiPos.Z_est.*10^-6;
                  %MultiPos.Yaw_est= MultiPos.Yaw_est.*10^-6;
                  %MultiPos.Pitch_est= MultiPos.Pitch_est.*10^-6;
                  %MultiPos.Roll_est= MultiPos.Roll_est.*10^-6;

                  %Only keep data from channel 6
                  %IDKeep = nan(length(NameMulti_Pos),1);

                  %for i3=1:length(NameMulti_Pos)
                      %St = NameMulti_Pos{i3};St=cellstr(St);St=St{1}(end-5:end);
                        %IDKeep(i3) = strcmp(St, '_6.tif');
                  %end

                  %if sum(IDKeep) ==0 %If only 5 bands => work with band 1
                        %for i3=1:length(NameMulti_Pos)
                            %St = NameMulti_Pos{i3};St=cellstr(St);St=St{1}(end-5:end);
                            %IDKeep(i3) = strcmp(St, '_1.tif');
                        %end
                  %end

                  %NameMulti_Pos = NameMulti_Pos(IDKeep==1); 
                  %MultiPos = MultiPos(IDKeep==1,:);
                  %clear IDKeep
                  %NameMulti_Pos = char(NameMulti_Pos);

                  %% Match multi lists
                  %[a,b] = ismember(NameMulti_Pos,ListMulti);b=unique(b);b(b==0)=[];
                  %NameMulti_Pos = NameMulti_Pos(a==1);
                  %MultiPos=MultiPos(a==1,:);
                  %ListMulti = ListMulti(b);
                  %TimeMulti = TimeMulti(b);
                  %LabelMulti =LabelMulti(b);


                   %TimeMulti_corr = TimeMulti;
                    %for i6=2:length(TimeMulti_corr)
                        %if TimeMulti_corr(i6)==TimeMulti_corr(i6-1)
                            %TimeMulti_corr(i6)=TimeMulti_corr(i6)+0.5./(24.*3600);
                        %end
                    %end

                  %Remove Non-Nan values
                  %ID = ~isnan(MultiPos.X_est);

                      %MultiPos = MultiPos(ID==1,:);
                      %NameMulti_Pos = NameMulti_Pos(ID==1);
                      %TimeMulti_corr=TimeMulti_corr(ID==1);

                  %clear ID

                  %% Match the data
                        %NameThermal = struct2cell(DirIm);NameThermal=NameThermal';NameThermal = NameThermal(:,1);
                        %X = nan(length(NameThermal),1);Y=X;Z=X;yaw = X; pitch=X; roll=X;
                     %TableThermal=table(NameThermal,X,Y,Z,yaw,pitch,roll);clear X Y pitch yaw Z roll
                     %for i4=1:length(TableThermal.X)
                         %try
                            %TableThermal.X(i4) = interp1(TimeMulti_corr, MultiPos.X_est, TTime_corr(i4)-TimeDiff/(24*3600)) + Bias(1);
                            %TableThermal.Y(i4) = interp1(TimeMulti_corr(:), MultiPos.Y_est, TTime_corr(i4)-TimeDiff/(24*3600))+ Bias(2);
                            %TableThermal.Z(i4) = interp1(TimeMulti_corr(:), MultiPos.Z_est, TTime_corr(i4)-TimeDiff/(24*3600))+ Bias(3);
                            %TableThermal.yaw(i4) = interp1(TimeMulti_corr, MultiPos.Yaw_est, TTime_corr(i4)-TimeDiff/(24*3600))+ Bias(4);
                            %TableThermal.pitch(i4) = interp1(TimeMulti_corr, MultiPos.Pitch_est, TTime_corr(i4)-TimeDiff/(24*3600))+ Bias(5);
                            %TableThermal.roll(i4) = interp1(TimeMulti_corr, MultiPos.Roll_est, TTime_corr(i4)-TimeDiff/(24*3600))+ Bias(6);
                         %end
                     %end

                     %Look for NaN-values, and replace them with mean value
                     %from above or below
                     %A = find(isnan(TableThermal.X));
                     %if ~isempty(A)
                         %for k1=1:length(A)  
                            %try
                            %TableThermal.X(A(k1)) = nanmean(TableThermal.X(max(1,A(k1)-5):min(A(k1)+5, length(TableThermal.X))));
                            %TableThermal.Y(A(k1)) = nanmean(TableThermal.Y(max(1,A(k1)-5):min(A(k1)+5, length(TableThermal.X))));
                            %TableThermal.Z(A(k1)) = nanmean(TableThermal.Z(max(1,A(k1)-5):min(A(k1)+5, length(TableThermal.X))));
                            %TableThermal.yaw(A(k1)) = nanmean(TableThermal.yaw(max(1,A(k1)-5):min(A(k1)+5, length(TableThermal.X))));
                            %TableThermal.pitch(A(k1)) = nanmean(TableThermal.pitch(max(1,A(k1)-5):min(A(k1)+5, length(TableThermal.X))));
                            %TableThermal.roll(A(k1)) = nanmean(TableThermal.roll(max(1,A(k1)-5):min(A(k1)+5, length(TableThermal.X))));

                        %end
                     %end


              %Handles.H=figure('Name','XY-pattern Thermal');  scatter(TableThermal.X,TableThermal.Y, 'o');
                %prompt = 'Correct XY-values? (based on multispectral data) (1/0)';
                %x = input(prompt);
                %close(Handles.H);
            %if x==1
                %writetable(TableThermal,[LocExp 'ThermalPos_FromMulti.txt']);
            %else
                %display('Issue with matching thermal and multispectral image position: Data not saved');
            %end


            %else
                %display('Issue with matching thermal and multispectral image position: metadata thermal not correct');
                     %end
            %end
         %catch
              %display('Issue with matching thermal and multispectral image position')
         %end
    %end

%if no data from multispectral flights are available, get XYZ from the GNSS of the thermal camera.
    if length(DirIm) == length(ThermalInf.Date)
            NameThermal = struct2cell(DirIm);NameThermal=NameThermal';NameThermal = NameThermal(:,1);
            
            TableThermal=table(NameThermal,ThermalInf.Latitude,ThermalInf.Longitude,ThermalInf.Altitude,ThermalInf.Yaw,ThermalInf.Pitch,ThermalInf.Roll, 'VariableNames', {'Names', 'lat','lon','alt', 'yaw','pitch','roll',});
                        
            Handles.H=figure('Name','XY-pattern Thermal');  scatter(TableThermal.lat,TableThermal.lon, 'o');
            prompt = 'Correct XY-values? (based on thermal GNSS sensor) (1/0)';
            x = input(prompt);
            close(Handles.H);
            writetable(TableThermal,[LocExp 'ThermalPos_FromThermal.txt']);
    else
       disp('Issue with GNSS position thermal camera: metadata thermal not correct');
    end



%% 

function [sharpness]=estimate_sharpness(G)

[Gx, Gy]=gradient(G);

S=sqrt(Gx.*Gx+Gy.*Gy);
sharpness=sum(sum(S))./(numel(Gx));

end

%%
function A= f_readIm(ImName, Opt_CSV)

  %Open data
   if Opt_CSV
     A=readmatrix(ImName);
   else
    A = imread(ImName);A=double(A); A=A.*0.04-273.15;
   end
end

% 
% function namemeta = ImportMeta(filename)
% %IMPORTFILE Import data from a text file
% %  NAMEMETA = IMPORTFILE(FILENAME) reads data from text file FILENAME
% %  for the default selection.  Returns the data as a table.
% %
% %  NAMEMETA = IMPORTFILE(FILE, DATALINES) reads data for the specified
% %  row interval(s) of text file FILENAME. Specify DATALINES as a
% %  positive scalar integer or a N-by-2 array of positive scalar integers
% %  for dis-contiguous row intervals.
% %
% %  Example:
% %  namemeta = importfile("E:\2021_drone\20210814_Gedinne\Gedinne1\thermal\testje\name_meta.csv", [2, Inf]);
% %
% %  See also READTABLE.
% %
% % Auto-generated by MATLAB on 25-Nov-2021 12:04:15
% 
% %% Input handling
% 
% % If dataLines is not specified, define defaults
% %if nargin < 2
%     dataLines = [2, Inf];
% %end
% 
% %% Setup the Import Options and import the data
% opts = delimitedTextImportOptions("NumVariables", 14);
% 
% % Specify range and delimiter
% opts.DataLines = dataLines;
% opts.Delimiter = ";";
% 
% % Specify column names and types
% opts.VariableNames = ["Date", "Time", "Latitude", "Longitude", "Pitch", "Roll", "Yaw", "Altitude", "Satellites", "TrackHeading", "MagneticHeading", "TempHousing", "TempFPA", "FFC_AfterFrame"];
% opts.VariableTypes = ["categorical", "datetime", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
% 
% % Specify file level properties
% opts.ExtraColumnsRule = "ignore";
% opts.EmptyLineRule = "read";
% 
% % Specify variable properties
% opts = setvaropts(opts, "Date", "EmptyFieldRule", "auto");
% opts = setvaropts(opts, "Time", "InputFormat", "HH:mm:ss");
% 
% % Import the data
% namemeta = readtable(filename, opts);
% 
% end


function namemeta = ImportMeta(filename)
%IMPORTFILE Import data from a text file
%  NAMEMETA = IMPORTFILE(FILENAME) reads data from text file FILENAME
%  for the default selection.  Returns the data as a table.
%
%  NAMEMETA = IMPORTFILE(FILE, DATALINES) reads data for the specified
%  row interval(s) of text file FILENAME. Specify DATALINES as a
%  positive scalar integer or a N-by-2 array of positive scalar integers
%  for dis-contiguous row intervals.
%
%  Example:
%  namemeta = importfile("E:\2021_drone\20210814_Gedinne\Gedinne1\thermal\testje\name_meta.csv", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 25-Nov-2021 12:04:15

%% Input handling

% If dataLines is not specified, define defaults
%if nargin < 2
    dataLines = [2, Inf];
%end

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 14);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ";";

% Specify column names and types
opts.VariableNames = ["Date", "Time", "Latitude", "Longitude", "Pitch", "Roll", "Yaw", "Altitude", "Satellites", "TrackHeading", "MagneticHeading", "TempHousing", "TempFPA", "FFC_AfterFrame"];
opts.VariableTypes = ["categorical", "datetime", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "Date", "EmptyFieldRule", "auto");
opts = setvaropts(opts, "Time", "InputFormat", "HH:mm:ss");

% Import the data
namemeta = readtable(filename, opts);

end

%function NameMulti_Pos =Imp_MultiNames(filename)
%IMPORTFILE Import data from a text file
%  GED1MULTIPOS = IMPORTFILE(FILENAME) reads data from text file
%  FILENAME for the default selection.  Returns the data as a cell array.
%
%  GED1MULTIPOS = IMPORTFILE(FILE, DATALINES) reads data for the
%  specified row interval(s) of text file FILENAME. Specify DATALINES as
%  a positive scalar integer or a N-by-2 array of positive scalar
%  integers for dis-contiguous row intervals.
%
%  Example:
%  Ged1MultiPos = importfile("E:\2021_drone\20210814_Gedinne\Agisoft\210814_Ged1_Multi_Pos.txt", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 25-Nov-2021 12:51:37

%% Input handling

% If dataLines is not specified, define defaults

    %dataLines = [2, Inf];


%% Setup the Import Options and import the data
%opts = delimitedTextImportOptions("NumVariables", 13);

% Specify range and delimiter
%opts.DataLines = dataLines;
%opts.Delimiter = ",";

% Specify column names and types
%opts.VariableNames = ["Label", "Var2", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13"];
%opts.SelectedVariableNames = "Label";
%opts.VariableTypes = ["string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string"];

% Specify file level properties
%opts.ExtraColumnsRule = "ignore";
%opts.EmptyLineRule = "read";

% Specify variable properties
%opts = setvaropts(opts, ["Label", "Var2", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13"], "WhitespaceRule", "preserve");
%opts = setvaropts(opts, ["Label", "Var2", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13"], "EmptyFieldRule", "auto");

% Import the data
%Ged1MultiPos = readtable(filename, opts);

%% Convert to output type
%NameMulti_Pos = table2cell(Ged1MultiPos);
% numIdx = cellfun(@(x) ~isnan(str2double(x)), Ged1MultiPos);
% NameMulti_Pos(numIdx) = cellfun(@(x) {str2double(x)}, Ged1MultiPos(numIdx));
%end

%function MultiPos = Import_Pos(filename)
%IMPORTFILE Import data from a text file
%  GED1MULTIPOS = IMPORTFILE(FILENAME) reads data from text file
%  FILENAME for the default selection.  Returns the data as a table.
%
%  GED1MULTIPOS = IMPORTFILE(FILE, DATALINES) reads data for the
%  specified row interval(s) of text file FILENAME. Specify DATALINES as
%  a positive scalar integer or a N-by-2 array of positive scalar
%  integers for dis-contiguous row intervals.
%
%  Example:
%  Ged1MultiPos = importfile("E:\2021_drone\20210814_Gedinne\Agisoft\210814_Ged1_Multi_Pos.txt", [1, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 25-Nov-2021 12:38:26

%% Input handling

% If dataLines is not specified, define defaults

    %dataLines = [2, Inf];


%% Setup the Import Options and import the data
%opts = delimitedTextImportOptions("NumVariables", 13);

% Specify range and delimiter
%opts.DataLines = dataLines;
%opts.Delimiter = ",";

% Specify column names and types
%opts.VariableNames = ["Var1", "Var2", "Var3", "Var4", "Var5", "Var6", "Var7", "X_est", "Y_est", "Z_est", "Yaw_est", "Pitch_est", "Roll_est"];
%opts.SelectedVariableNames = ["X_est", "Y_est", "Z_est", "Yaw_est", "Pitch_est", "Roll_est"];
%opts.VariableTypes = ["string", "string", "string", "string", "string", "string", "string", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
%opts.ExtraColumnsRule = "ignore";
%opts.EmptyLineRule = "read";

% Specify variable properties
%opts = setvaropts(opts, ["Var1", "Var2", "Var3", "Var4", "Var5", "Var6", "Var7"], "WhitespaceRule", "preserve");
%opts = setvaropts(opts, ["Var1", "Var2", "Var3", "Var4", "Var5", "Var6", "Var7"], "EmptyFieldRule", "auto");
%opts = setvaropts(opts, ["X_est", "Y_est", "Z_est", "Yaw_est", "Pitch_est", "Roll_est"], "TrimNonNumeric", true);
%opts = setvaropts(opts, ["X_est", "Y_est", "Z_est", "Yaw_est", "Pitch_est", "Roll_est"], "DecimalSeparator", ",");
%opts = setvaropts(opts, ["X_est", "Y_est", "Z_est", "Yaw_est", "Pitch_est", "Roll_est"], "ThousandsSeparator", ".");

% Import the data
%MultiPos = readtable(filename, opts);

%end
