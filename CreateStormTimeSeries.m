%% ------------------------CreateStormTimeSeries---------------------------%
% Purpose: This script creates a time series of storms
% from hindacasts of water levels and wave conditions.
%
% References:
% ...[1] Wahl, T., Plant, N. G., & Long, J. W. (2016). Probabilistic 
% assessment of erosion and flooding risk in the northern Gulf of Mexico. 
% Journal of Geophysical Research: Oceans, 121(5), 3029-3043.
%
%       IRBReeves    8/4/23
%       WBP Franklin 8/18/24 (Added functionality to export generated TS)
%
% Variable naming rules:
%       r - real array
%       n - integer value
%       f - float value
%       b - boolean
%       s - string
%       dt - datetime
%
% Inputs:
%       input_filename       - filename of hourly wave & water level records
%       fBeta                - beach slope 
%       fBermEl              - berm elevation (m NAVD88)
%       bPlot                - boolean for plotting
%
%% -----------------------------------------------------------------------%

clear

% Define Inputs
%input_filename = '/Users/reevesi/Desktop/WHOI-USGS/Data/Storms/NCoreBanks_WL_Hs_Tp.mat';  % '/Users/reevesi/Desktop/WHOI-USGS/Data/Storms/NCBdata_nans.mat';
input_filename = 'C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 3\\Ian_MSSM_RENCI\\MSSM_RENCI\\MSSM_RENCI\\NcoreBanks_WL_Hs_Tp.mat'
fBeta = 0.039;
fBermEl = 1.78;  % NCB: [Avg] 2.03, [25%] 1.78, [%10] 1.56, [%5] 1.45
fTWL_thresh = 0.39; % MHW NAVD88
site = 4;  % Reference number for TWL/Wave timeseries station
storm_step = 25; % Number of storm storm iterations in a year
model_step = 50; % Number of iterations in model year
bPlot = true;
MHW = 0.26; % MHW for Core Banks area (from USGS Report, Morton and Miller; 2005)

% Run Storm Time Series Function
[stStorms, stHWH, storm_time_series, stIter_Storm_Prob] = StormTimeSeries(input_filename, fBeta, fBermEl, fTWL_thresh, site, storm_step, ...
    model_step, bPlot,MHW);


%% ------------------------------------------------------------------------%
function [stStorms, stHWH, storm_time_series, stIter_Storm_Prob] = StormTimeSeries(input_filename, fBeta, fBermEl, fTWL_thresh, site, storm_step, ...
    model_step, bPlot,MHW)

% thresholds for identifying storms
nMinDuration = 8; % hr, minimum duration of a storm 

nStartYear = 1979; % just used for plotting
    
%% -----------------------------------------------------------------------%

% load USGS data
stObs = load_data_USGS_NCB(input_filename, site); 

% calculate the TWL and R2% from observational record
stObs = calculate_TWL(fBeta, stObs);

% extract sea-storm events from the observational record
stStorms = extract_sea_storms_from_obs(fBermEl, nMinDuration, stObs);

% extract high water hours (HWHs) from obervational record
stHWH = extract_HWH_from_obs(fTWL_thresh,model_step,stObs);

% create time series from observed storms
storm_time_series = create_time_series(stObs, stStorms, storm_step, model_step);

% extract probability of storm occurance for each iteration of a year
stIter_Storm_Prob = extract_iteration_storm_probabilities(storm_time_series, model_step, storm_step);

export_storms_for_Cascade = export_storm_TS(MHW,stStorms)

% Print number of storms
stStorms.nStorms

%% -----------------------------------------------------------------------%

    function stObs = load_data_USGS_NCB(input_filename, site)
    
        % save as structure
        stObs = struct();
        
        % USGS simulated wave and water level hindcast for North Core Banks, NC, 1980-2020    
        
        rData = load(input_filename);
        stObs.dtH   = datetime(rData.TT, 'ConvertFrom', 'datenum');
        stObs.rHs   = rData.HS{site};
        stObs.rTp   = rData.TP{site};  % 11 = TPD, 12 = TP
        stObs.rWavD = rData.WDIR{site};
        stObs.rSL   = rData.ZZ{site};
        stObs.dtSL  = stObs.dtH;

    end

    function stObs = calculate_TWL(fBeta, stObs)

        % calculate R2% and add to SL to get the TWL (currently only corrected data)     
    
        % KA: is this Stockdon 2006 broken down into components?
        stObs.rL0   = (9.8 * stObs.rTp.^2) / (2 * pi); % wavelength  
        rSetup = 0.35 * fBeta * sqrt(stObs.rHs .* stObs.rL0); 
        rSin   = 0.75 * fBeta * sqrt(stObs.rHs .* stObs.rL0);  % incident band swash
        rSig   = 0.06 * sqrt(stObs.rHs .* stObs.rL0) ;         % infragravity band swash
        rSwash = sqrt((rSin.^2) + (rSig.^2));      % total swash
        stObs.rR2    = 1.1 * (rSetup + (rSwash/2));      % R2%
    
        stObs.rTWL  = stObs.rSL + stObs.rR2;
        stObs.rRlow = (stObs.rTWL - (rSwash/2));
    
        if bPlot == 1       
            figure
            % SL
            subplot(3,2,1)
            plot(stObs.dtH, stObs.rSL)
            ylabel('Sea level [mNAVD88]')
    
            % Hs
            subplot(3,2,2)
            plot(stObs.dtH, stObs.rHs)
            ylabel('Hs [m]')
            ylim([0 10])
    
            % Tp
            subplot(3,2,3)
            plot(stObs.dtH, stObs.rTp)
            ylabel('Tp [s]')
            ylim([0 20])
    
            % WavD
            subplot(3,2,4)
            plot(stObs.dtH, stObs.rWavD)
            ylabel('Wave Direction [degree]')
            
            % R2
            subplot(3,2,5)
            plot(stObs.dtH, stObs.rR2)
            ylabel('R2/% [m]')
    
            % TWL
            subplot(3,2,6)
            plot(stObs.dtH, stObs.rTWL)
            ylabel('TWL [m NAVD88]')
        end
    end

    function stStorms = extract_sea_storms_from_obs(fBermEl,nMinDuration,stObs)

        %%% First, define TWL threshold exceedence. Then, how much the different 
        %%% variables contributed to those events and if there is a dominant driver 
        %%% that can be used for the event selection

        % for each year, find when the TWL exceeds an erosion threshold

        % find the annual average TWL from all threshold exceedances from a 
        % given year calculate annual averages of MSL (here 35-day mean), 
        %  and R2% during the TWL exceedances

        stObs.dtYears = year(stObs.dtH);  % array of year for each hourly measurement

        nYrs = max(stObs.dtYears) - min(stObs.dtYears) + 1;
        [rHs_over_yearly, rR2_over_yearly, rTWL_over_yearly] = deal(NaN(nYrs,1));

        year_start = stObs.dtYears(1);
        year_stop = stObs.dtYears(end);

        for iYear = year_start : year_stop

            idx = find(ismember(stObs.dtYears, iYear) == 1);
            iStart = idx(1);
            iStop = idx(end);

            rHH = stObs.rHs(iStart:iStop);      
            rRR = stObs.rR2(iStart:iStop); 
            rTT = stObs.rTWL(iStart:iStop);     
            rTWL_over_yearly(iYear - year_start + 1) = mean(rTT(rTT > fBermEl), "omitnan");
            rHs_over_yearly(iYear - year_start + 1)  = mean(rHH(rTT > fBermEl), "omitnan");
            rR2_over_yearly(iYear - year_start + 1)  = mean(rRR(rTT > fBermEl), "omitnan");

        end

        % identify Hs threshold to qualify as a storm event, round nearest 0.05 m
        % nHs_min = min(rHs_over_yearly);
        nHs_min = prctile(rHs_over_yearly, 5);
        nHs_threshold = floor(nHs_min / 0.05) * 0.05;     



        % Alternative method: Define Hs from lower 2 sigma of all TWL exceedances
        % nHs_min = mean(rHH(rTT > rBermEl)) - 2 * std(rHH(rTT > rBermEl));
        % nHs_threshold = floor(nHs_min / 0.05) * 0.05

        % visual check of threshold and drivers (this is hard coded to 1980, will 
        % need to be changed for other locations)
        if bPlot == 1
            figure
            subplot(2,1,1)
            plot(nStartYear : 1 : nStartYear+nYrs-1, [rTWL_over_yearly, rR2_over_yearly], '-o')
            ylabel('Contribution to TWL')
            legend('TWL', 'R2/%')
            subplot(2,1,2)
            plot(nStartYear : 1 : nStartYear+nYrs-1, rHs_over_yearly, '-o')
            hold on
            refline(0, nHs_threshold) % * Requires Statistics & Machine Learning toolbox *
            ylabel('Hs')
            title(sprintf('Berm Elev = %d m NAVD88', fBermEl))
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % find storms (translated from Ian Reeves); use corrected data only
        [iStormStart, iStormStop, dtStormStart, dtStormStop, cStormHs, ...
            cStormDur, cStormTWL, cStormRlow, cStormSL, cStormTp, ...
            cStormWavD, cStormNegSurgeDT, dtYear] = deal(cell(0));

        % VCR storms are all significantly under these max thresholds
        dtH_yr = year(stObs.dtH);
        t = 1 ;

        while t <= length(stObs.dtH)

            if stObs.rHs(t) >= nHs_threshold
                stormStart = t;
                dur = 1;
                t = t + 1;

                % if Hs drops below Hs_threshold for only 24 hrs or less, 
                % exceedence is assumed part of same weather system 
                % (Wahl et al., 2016; Li et al., 2014)

                try
                    while t < (length(stObs.rHs) - 24) && sum(stObs.rHs(t:t+24) > nHs_threshold) > 0 
                        if stObs.rHs(t) > nHs_threshold
                            dur = dur + 1;
                            t = t + 1;
                        else
                            t = t + 1;
                        end
                    end
                catch
                    keyboard
                end

                % minimum of an 8 hr storm (Magliocca et al., 2011)
                if dur > nMinDuration
                    stormStop = t;
                    iStormStart{end+1} = stormStart;
                    iStormStop{end+1} = stormStop;
                    dtStormStart{end+1} = datenum(stObs.dtH(stormStart));
                    dtStormStop{end+1} = datenum(stObs.dtH(stormStop));
                    cStormDur{end+1} = dur;
                    cStormTWL{end+1} = max(stObs.rTWL(stormStart:stormStop));
                    cStormRlow{end+1} = max(stObs.rRlow(stormStart:stormStop));
                    cStormSL{end+1} = max(stObs.rSL(stormStart:stormStop)); % Sea Level (NTR + AT)

                    % need to find the max Hs and simultaneous
                    % (not max Tp and WavD)
                    [cStormHs{end+1}, iHs] = max(stObs.rHs(stormStart:stormStop));
                    tmpTp = stObs.rTp(stormStart:stormStop);
                    tmpWavD = stObs.rWavD(stormStart:stormStop);
                    cStormTp{end+1} = tmpTp(iHs);
                    cStormWavD{end+1} = tmpWavD(iHs);
                    dtYear{end+1} = dtH_yr(stormStart);
                end

                t = t + 1;  
            else

                t = t + 1;

            end
        end

        % create structure of storm parameters
        stStorms = struct();
        stStorms.rStart_dt = cell2mat(dtStormStart)';
        stStorms.rStop_dt = cell2mat(dtStormStop)';
        stStorms.rStart = datetime(stStorms.rStart_dt, 'ConvertFrom', 'datenum');
        stStorms.rStop = datetime(stStorms.rStop_dt, 'ConvertFrom', 'datenum');
        stStorms.rHs = cell2mat(cStormHs)';
        stStorms.rDur = cell2mat(cStormDur)';
        stStorms.rTWL = cell2mat(cStormTWL)';
        stStorms.rSL = cell2mat(cStormSL)';
        stStorms.rTp = cell2mat(cStormTp)';
        stStorms.rRlow = cell2mat(cStormRlow)';
        stStorms.rYear = cell2mat(dtYear)';
        stStorms.nHs_threshold = nHs_threshold;
        stStorms.nStorms = length(stStorms.rTWL);
        stStorms.rWeek = week(stStorms.rStart);
        stStorms.rDay = day(stStorms.rStart, 'dayofyear');

        if bPlot == 1
            % Plot storm TWL, Hs, Dur, Tp, NTR, AT histogram
            figure
            subplot(2,2,1)
            histogram(stStorms.rTWL, 50)
            ylabel('Storm TWL [m NAVD88]')

            subplot(2,2,2)
            histogram(stStorms.rHs, 50)
            ylabel('Storm Hs [m]')
            %title('1980 - 2014')

            subplot(2,2,3)
            histogram(stStorms.rDur, 50)
            ylabel('Storm Dur [hrs]')

            subplot(2,2,4)
            histogram(stStorms.rTp, 50)
            ylabel('Storm Tp [s]')
        end

    end

    function storm_time_series = create_time_series(stObs, stStorms, storm_step, model_step)

        num_storms = length(stStorms.rTWL);
        start_year = min(year(stObs.dtH));
        [iter, rhigh, rlow, dur] = deal(cell(0));

        for s = 1 : num_storms
            step = floor(stStorms.rDay(s) * (storm_step / 365));     
            iteration = (stStorms.rYear(s) - start_year) * storm_step + step;

            if s > 1 && iteration == prev_iter
                if stStorms.rTWL(s) > prev_twl
                    % Store new in place of prev
                    iter{end} = iteration;
                    rhigh{end} = stStorms.rTWL(s);
                    rlow{end} = stStorms.rRlow(s);
                    dur{end} = min(stStorms.rDur(s), floor(365 / storm_step * 24)); % Maximum storm duration is length [hrs] of storm iteration

                    % Store previous
                    prev_iter = iteration;
                    prev_twl = stStorms.rTWL(s);
                end  % Don't store storm in same storm_step if of lesser TWL
            else
                %store
                iter{end + 1} = iteration;
                rhigh{end + 1} = stStorms.rTWL(s);
                rlow{end + 1} = stStorms.rRlow(s);
                dur{end + 1} = stStorms.rDur(s);

                % Store previous
                prev_iter = iteration;
                prev_twl = stStorms.rTWL(s);
            end
        end

        storm_time_series = NaN(length(rhigh), 4);
        storm_time_series(:, 1) = cell2mat(iter) .* (model_step / storm_step);
        storm_time_series(:, 2) = cell2mat(rhigh);
        storm_time_series(:, 3) = cell2mat(rlow);
        storm_time_series(:, 4) = cell2mat(dur);

    end

    function stHWH = extract_HWH_from_obs(fTWL_thresh,model_step,stObs)

        % Extract High Water Hours (HWH, hours in which TWL exceeds reference elevation) from observations

        stObs.dtDaysDif = daysdif(stObs.dtH, stObs.dtH(1)) * -1; % Number of days from start of time series  TODO: Need to insure start of time series is a multiple 50th of year
        % * Requires Statistics & Machine Learning toolbox *

        model_step_length = (365 * 24) / model_step;

        step_num = floor(length(stObs.dtDaysDif) / model_step_length);

        stHWH = NaN(step_num,3); % Initialize

        step = 1;

        for iStep = 1 : model_step_length : ceil(length(stObs.dtDaysDif)) - model_step_length

            iStart = floor(iStep);
            iStop = iStep + floor(model_step_length);

            rTT = stObs.rTWL(iStart:iStop); % All TWLs for this time step
            rTT_over = rTT(rTT > fTWL_thresh); % All TWLs for this time step in exceedance of threshold

            rTWL_over_dur = length(rTT_over); % Number of hours TWL exceeded threshold for this time step

            if rTWL_over_dur > 0
                rTWL_over_avg = mean(rTT_over, "omitnan"); % Average TWL for this time step for hours in which TWL exceeded threshold
                rTWL_over_max = max(rTT_over); % Maximum TWL for this time step for hours in which TWL exceeded threshold
            else
                rTWL_over_avg = 0;
                rTWL_over_max = 0;
            end

            stHWH(step, :) = [rTWL_over_dur, rTWL_over_max, rTWL_over_avg];

            step = step + 1;

        end

    end

    function stIter_Storm_Prob = extract_iteration_storm_probabilities(storm_time_series, model_step, storm_step)

        storm_iter = rem(storm_time_series(:,1), model_step);
        storm_iter_bins = (0:storm_step) * (model_step / storm_step);  % [days] Bin edges for when observed storms occured
        
        h = histogram(storm_iter, storm_iter_bins);  % Bin observed storms by iteration in which they occured

        bin_val = h.Values;  % Get bin values

        obs_dur = stStorms.rYear(end) - stStorms.rYear(1);  % [years] Duration of observations

        stIter_Storm_Prob = bin_val ./ obs_dur;  % Get probability of storm occurance for each iteration of the year

    end

    function save_TS = export_storm_TS(MHW,stStorms)
        
        % Seperate out the Rhigh and Rlow values
        Rhigh_m = stStorms.rTWL
        Rlow_m = stStorms.rRlow

        % Convert the Rhigh and Rlow values from meters to decameters
        Rhigh_dm = (Rhigh_m-MHW)/10
        Rlow_dm = (Rlow_m-MHW)/10
        
        % Assign model year
        Raw_Year_Data = stStorms.rYear
        Export_Year = Raw_Year_Data - stStorms.rYear(1)
        
        % Period
        Export_Period = stStorms.rTp

        % Storm Duration
        Raw_Duration = stStorms.rDur
        Export_Duration = round(Raw_Duration/2)

        % Combine all data together 
        combined_Matrix = [Export_Year,Rhigh_dm,Rlow_dm,Export_Period,Export_Duration]
        writematrix(combined_Matrix,'Core_Banks_Historic_Storm_Record.csv')
        save_TS = combined_Matrix

    end

end