function [ records header ] = ncs_wrapper( ncs_path, sample_range, channels  )
%NCS_WRAPPER Wrapper for Neuralynx NCS format
%   Encapsulates the reading of Neuralunx NCS files. The name of the
%   channel recorded in the file is assumed to be presented in the last
%   numeric value of the file name, preceded by -, e.g. Cage1-1.ncs to
%   indicate channel 1 from cage 1.
%   Parameters:
%       ncs_path: Path to folder containing the NCS files. If NCS file is
%                 given, the path to the file is used.
%       sample_range 1x2 (optional) : Range of samples to be extracted.
%                                     Default: All samples.
%       channels 1xN (optional) : Channels to be extracted
%                                 Default: All channels
%  Output:
%       records: ChannelsxSamples matrix
%       header: Struct with fields:
%                   time_created : File creation time
%                   time_closed : File closing time
%                   sampling_frequency : EEG sampling frequency
%  Required files:
%  

    [path,name,ext] = fileparts(ncs_path);
    
    if isempty(name)
       ncs_path = path; 
    end

    SAMPLES_PER_RECORD = 512;
    
    extraction_mode = 1;
    extraction_vector = 1;
    if nargin > 1 && ~isempty(sample_range)
        extraction_mode = 2;
        extraction_vector = [1 max(floor(sample_range(2)),1)];       
    end
    
    records = [];   
    header = struct;    
    files = dir(fullfile(ncs_path, '*.ncs'));
     
    if nargin < 3
        channels = 1:size(files,1);
    end
    header.records = 0;
    header.duration = 0;
    header.label = cell(1,size(files,1));
    for c = 1:size(files,1)

        full_file = strcat(ncs_path ,filesep, files(c).name);
        %Extract selected range from the current NCS file
        if isunix 
            [ timestamps samples, ncs_header]  = Nlx2MatCSC_v3(full_file, [1,0,0,0,1],1,extraction_mode,extraction_vector);
        else
            [ timestamps samples, ncs_header]  = Nlx2MatCSC(full_file, [1,0,0,0,1],1,extraction_mode,extraction_vector);       
        end
        ADBitVolts = regexp(ncs_header(17,:),'(?<=(-ADBitVolts\s))[\s.0-9]*$','match');    
        ADBitVolts = str2double(ADBitVolts{1});
        samplingFrequency = regexp(ncs_header(15,:),'(?<=(-SamplingFrequency\s))[\s.0-9]*$','match');   
        samplingFrequency = str2double(samplingFrequency{1});          
        n_samples = size(samples,2)*SAMPLES_PER_RECORD;
           
        %Channel number is assumed to be presented as the last numeric
        %value of the file name.
        channel_number = str2double(regexp(files(c).name, '(?<=[-])[\d]*(?=(.ncs))', 'match'));
        
        %Skip unselected channels        
        if ~ismember(channel_number,channels)
           continue 
        end            
        ADBitVolts = regexp(ncs_header(17,:),'(?<=(-ADBitVolts\s))[\s.0-9]*$','match');    
        ADBitVolts = str2double(ADBitVolts{1});    
        channel_index = find( channels == channel_number);
        if isempty(records)
           records = zeros(size(channels,2), n_samples ); 
           header.time_created = regexp(ncs_header(8,:),'(?<=(-TimeCreated\s))[\s0-9:\\/]*$','match');
           header.time_closed = regexp(ncs_header(9,:),'(?<=(-TimeClosed\s))[\s0-9:\\/]*$','match');    
           frequency = regexp(ncs_header(15,:),'(?<=(-SamplingFrequency\s))[\s0-9]*$','match');
           header.frequency = str2double(frequency{1});   
           header.ADBitVolts = ADBitVolts;
        end
        label = regexp(ncs_header(15,:),'(?<=(-AcqEntName\s))[\s0-9:A-Za-z\\/\-]*$','match');
        header.label{channel_index} = char(label{1});     
        records(channel_index,1:n_samples) = reshape(samples,1,n_samples);
        records(channel_index,:) = ADBitVolts*records(channel_index,:);      
        header.records = max(header.records, n_samples);
        header.duration = max(header.duration, n_samples/header.frequency);
    end
end


