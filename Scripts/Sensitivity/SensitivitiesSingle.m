% Script tu run the sensitivity_analyses.m function using our data
% extracted from the script "SenseFileConv.R" to compute the sensitivity
% indexses for the single experiment inference results.

files4 = [string('Calibration_4'),string('Calibration_5'),...
    string('Calibration_6'),string('DynStim_1'), string('DynStim_2'), ...
    string('DynStim_3'), string('DynStim_8'), string('DynStim_9'),...
    string('DynStim_11'), string('DynStim_14')];

% Sensitivity index calculation for RFP results
for fn=1:length(files4)
    
    % Extract for each iteration the corresponding parameter draws matrix
    s1 = strcat('draws_', files4(fn), '.mat');
    post=load(s1{1});
    post = cell2mat(struct2cell(post));
    [nSamples,nPars]=size(post);
    
    % Define number of bins (square root of the number of samples)
    nBins=90;
    
    % Extract for each iteration the corresponding predictions matrix for
    % the samples extracted
    s2 = strcat('outpR_', files4(fn), '_RFP.mat');
    outp=load(s2{1});
    outp = cell2mat(struct2cell(outp));
    [~,nPoints] = size(outp);
    
    % Compute sensitivity indexces using the Olivia Eriksson function
    [S_i, all_par_vars]=sensitivity_analyses(post,outp, nBins);
    
    % Save the results in a csv file
    cr1 = strcat('RFPsense_',files4(fn), '.csv' );
    csvwrite(cr1{1}, S_i);
end

% Sensitivity index calculation for GFP results
for fn=1:length(files4)
    
    % Extract for each iteration the corresponding parameter draws matrix
    s1 = strcat('draws_', files4(fn), '.mat');
    post=load(s1{1});
    post = cell2mat(struct2cell(post));
    [nSamples,nPars]=size(post);

    % Define number of bins (square root of the number of samples)
    nBins=90;
    
    % Extract for each iteration the corresponding predictions matrix for
    % the samples extracted
    s2 = strcat('outpG_', files4(fn), '_GFP.mat');
    outp=load(s2{1});
    outp = cell2mat(struct2cell(outp));
    [~,nPoints] = size(outp);
    
    % Compute sensitivity indexces using the Olivia Eriksson function
    [S_i, all_par_vars]=sensitivity_analyses(post,outp, nBins);
    
    % Save the results in a csv file
    cr2 = strcat('GFPsense_',files4(fn), '.csv' );
    csvwrite(cr2{1}, S_i);
end

% Extraction of theSensitivity index vectors for all experiments taking
% only the parameters of interest; k_aTc (rK and gK) and theta_L (rTL, gTL)
rK = NaN(10,271);
rTL = NaN(10,271);

gK = NaN(10,271);
gTL = NaN(10,271);

for fn=1:length(files4)
    
    % Load the results for both TetR and LacI
    s1 = strcat('RFPsense_', files4(fn), '.csv');
    
    post=load(s1{1});
    
    s2 = strcat('GFPsense_', files4(fn), '.csv');
    
    post2=load(s2{1});
    
    % introduce corresponding vectors into the empty matrices where NaN
    % values are predefined to account for the different length of the
    % different experimental time points. 
    f = size(post);
    for tp=1:f(2)
        
       rK(fn,tp) = post(2,tp);
       rTL(fn,tp) = post(11,tp);
       
       gK(fn,tp) = post2(2,tp);
       gTL(fn,tp) = post2(11,tp);
        
    end
    
    
end


% Save results in CSV files
csvwrite('Sense_RFP_KaTc_.csv', rK');
csvwrite('Sense_RFP_ThetaL_.csv', rTL');
csvwrite('Sense_GFP_KaTc_.csv', gK');
csvwrite('Sense_GFP_ThetaL_.csv', gTL');






