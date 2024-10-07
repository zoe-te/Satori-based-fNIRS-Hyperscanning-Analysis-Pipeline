%% Analysis for fNIRS Hyperscanning (Correlation coefficient, Fisher Z Transformation, t-test)

% Define the number of participants (adjust based on your data)
num_participants = 2;

% Define the number of conditions (adjust based on your data)
num_conditions = 2;

% Define the base file path (adjust to your file path)
base_path = 'C:\Users\azt19\OneDrive\Dokumente\MATLAB\SatoriBasedAnalysis\';

% Total number of rows per condition (each condition must have the same number of rows)
rows_per_condition = 33;  % Adjust this based on your data structure

% Initialize cell arrays to store Fisher Z-scores and correlations for each condition
z_scores = cell(num_conditions, 1);
correlations = cell(num_conditions, 1);

%% Main loop
% Loop through each condition
for condition = 1:num_conditions
    % Calculate the row range dynamically for each condition
    start_row = (condition - 1) * rows_per_condition + 1;
    end_row = condition * rows_per_condition;
    cond_name = sprintf('cond%d', condition);
    
    % Initialize an array to store Z-scores and correlations for this condition
    z_score_table = [];
    correlation_table_all = [];  % Store correlation results for all pairs in this condition
    
    % Loop through each participant
    for participant = 1:num_participants
        % Create the filename for the current participant (e.g., 'ERAdata_001.csv')
        filename = sprintf('ERAdata_%03d.csv', participant);
        file_path = fullfile(base_path, filename);
        
        % Load the data from the CSV file
        data = readmatrix(file_path);
        
        % Extract the relevant rows for the current condition
        Cond = data(start_row:end_row, :);
        
        % Create a new file name for the output (e.g., 'ERAdata_001_cond1.csv')
        output_filename = sprintf('ERAdata_%03d_%s.csv', participant, cond_name);
        output_file_path = fullfile(base_path, output_filename);
        
        % Write the extracted rows to a new CSV file
        writematrix(Cond, output_file_path);
    end

    % Handle the long format table and correlation for the pairs of participants
    for participant = 1:2:num_participants
        % Load the data for the first participant in the pair (e.g., 'ERAdata_001_cond1.csv')
        filename1 = sprintf('ERAdata_%03d_%s.csv', participant, cond_name);
        file_path1 = fullfile(base_path, filename1);
        data1 = readmatrix(file_path1);
        
        % Load the data for the second participant in the pair (e.g., 'ERAdata_002_cond1.csv')
        filename2 = sprintf('ERAdata_%03d_%s.csv', participant + 1, cond_name);
        file_path2 = fullfile(base_path, filename2);
        data2 = readmatrix(file_path2);
        
        % Get the number of channels and timepoints for both participants
        [num_channels1, num_timepoints1] = size(data1);
        [num_channels2, num_timepoints2] = size(data2);
        
        % Create the arrays for participant, channel, and timepoint for both participants
        participant_column1 = participant * ones(num_channels1 * num_timepoints1, 1);
        channel_column1 = repelem((1:num_channels1)', num_timepoints1);
        timepoint_column1 = repmat((1:num_timepoints1)', num_channels1, 1);
        value_column1 = reshape(data1', [], 1);

        participant_column2 = (participant + 1) * ones(num_channels2 * num_timepoints2, 1);
        channel_column2 = repelem((1:num_channels2)', num_timepoints2);
        timepoint_column2 = repmat((1:num_timepoints2)', num_channels2, 1);
        value_column2 = reshape(data2', [], 1);
        
        % Combine the data for both participants into one long format table
        participant_column = [participant_column1; participant_column2];
        channel_column = [channel_column1; channel_column2];
        timepoint_column = [timepoint_column1; timepoint_column2];
        value_column = [value_column1; value_column2];
        
        long_format_table = table(participant_column, channel_column, timepoint_column, value_column, ...
            'VariableNames', {'participant', 'channel', 'timepoint', 'value'});
        
        % Save the long format data to a new CSV file
        output_filename = sprintf('ERAdata_%03d_%03d_%s_long_format.csv', participant, participant + 1, cond_name);
        output_file_path = fullfile(base_path, output_filename);
        writetable(long_format_table, output_file_path);
        
        % Calculate the Pearson correlation coefficient for each channel
        unique_channels = unique(long_format_table.channel);
        correlation_coeffs = zeros(length(unique_channels), 1);

        % Loop over each channel and calculate correlation coefficients
        for i = 1:length(unique_channels)
            current_channel = unique_channels(i);
            channel_data = long_format_table(long_format_table.channel == current_channel, :);
            
            values_participant_1 = channel_data.value(channel_data.participant == participant);
            values_participant_2 = channel_data.value(channel_data.participant == participant + 1);
            
            correlation_coeffs(i) = corr(values_participant_1, values_participant_2);
        end
        
       % Store the correlation coefficients for the current condition
        correlation_table_all = [correlation_table_all; unique_channels, correlation_coeffs];

        % Perform Fisher Z-transformation
        z_values = 0.5 * log((1 + correlation_coeffs) ./ (1 - correlation_coeffs));
        z_score_table = [z_score_table; unique_channels, z_values];
    end
    
    % Save the Z-scores to a new table
    z_table = array2table(z_score_table, 'VariableNames', {'channel', 'z_score'});
    output_filename = sprintf('ERAdata_%s_z_scores.csv', cond_name);
    output_file_path = fullfile(base_path, output_filename);
    writetable(z_table, output_file_path);
    
    % Store the Z-scores for later analysis & store correlation coefficients
    z_scores{condition} = z_score_table(:, 2);
    correlations{condition} = correlation_table_all(:, 2); 
end

%% Visualization: Scatter Plot of Z-scores across conditions
if num_conditions > 1
    figure;
    scatter(z_scores{1}, z_scores{2}, 'filled');
    xlabel('Z-scores Condition 1');
    ylabel('Z-scores Condition 2');
    title('Scatter Plot of Z-scores between Condition 1 and Condition 2');
    grid on;
end

%% Statistical Analysis: Paired t-test on Z-scores
if num_conditions > 1
    [h, p_value, ci, stats] = ttest(z_scores{1}, z_scores{2});
    
    % Display results of the t-test
    fprintf('Paired t-test Results:\n');
    fprintf('p-value: %.4f\n', p_value);
    fprintf('t-statistic: %.4f\n', stats.tstat);
    fprintf('Degrees of freedom: %d\n', stats.df);
    fprintf('Confidence interval: [%.4f, %.4f]\n', ci(1), ci(2));

    % Save the t-test results to a file
    ttest_results = table(p_value, stats.tstat, stats.df, ci(1), ci(2), ...
        'VariableNames', {'p_value', 't_stat', 'df', 'ci_lower', 'ci_upper'});
    output_filename = 't_test_results.csv';
    output_file_path = fullfile(base_path, output_filename);
    writetable(ttest_results, output_file_path);
end