function [ records header ] = ncs_wrapper( ncs_path,varargin  )
%NCS_WRAPPER Wrapper for Neuralynx NCS format
%   Encapsulates the reading of Neuralynx NCS files. The name of the
%   channel recorded in the file is assumed to be presented in the last
%   numeric value of the file name, preceded by -, e.g. 'Cage1-1.ncs' to
%   indicate channel 1 from cage 1. Multiple series for each channel are
%   assumed to be denoted by _<integers> postfix, e.g. 'Cage1-1_0001.ncs'.
%   
%   Parameters:
%       ncs_path: 
%               Path to folder containing the NCS files. If NCS file is
%               given, the path to the file is used.
%       sample_range 1x2 (optional): 
%               Range of samples to be extracted.
%               Default: All samples.
%       channels 1xN (optional): 
%               Channels to be extracted
%               Default: All channels
%       callback function handle (optional): 
%               Callback function to be called after header extraction and
%               for each channel. Function must accept three parameters:
%               sample records (1xS), header struct and the index of the 
%               channel.
%               Default: None
%  Output:
%       records: ChannelsxSamples matrix
%       header: Struct with fields:
%                   time_created : File creation time (String)
%                   time_closed : File closing time (String)
%                   frequency : EEG sampling frequency for each channel 1xC
%                               matrix
%                   units : Channel units 1xC cell array.
%                   label : Channel label 1xC cell array.
%  Required files: Nlx2MatCSC_v3, Nlx2MatCSC
%  Author: rciszek



    is_function = @(f) isa(f, 'function_handle');

    parser = inputParser;
    addRequired(parser,'ncs_path',@ischar);
    addParameter(parser,'sample_range',[],@isnumeric);
    addParameter(parser,'channels',[],@isnumeric);
    addParameter(parser,'callback',false,is_function);

    parse(parser,ncs_path,varargin{:});

    ncs_path = parser.Results.ncs_path; 
    sample_range = parser.Results.sample_range; 
    channels = parser.Results.channels; 
    callback = parser.Results.callback;     

    [path,~,~] = fileparts(ncs_path);
   
    SAMPLES_PER_RECORD = 512;
    SAMPLE_RESOLUTION = 1000000;
    
    extraction_mode = 1;
    extraction_vector = 1;
    if ~isempty(sample_range)
        extraction_mode = 2;
        extraction_vector = [1 max(floor(sample_range(2)),1)];       
    end
    
    records = [];   
    header = struct;    
    files = dir(fullfile(ncs_path, '*.ncs'));
     
    header.records = 0;
    header.duration = 0;
    header.label = cell(1,size(files,1));
        
    file_map = containers.Map('KeyType','double','ValueType','any');

    %Map files for each channel
    for c = 1:size(files,1)   
        channel_number = str2double(regexp(files(c).name, '(?<=[-])[\d]*', 'match')); 
        if file_map.isKey(channel_number)
            channel_files = file_map(channel_number);
            channel_files = horzcat( channel_files, {files(c).name});
            file_map(channel_number) = channel_files;
        else
            file_map(channel_number) = {files(c).name};
        end
    end
    
    unique_channels = file_map.keys();   
    unique_channels = sort(cell2mat(unique_channels));
    
    if isempty(channels)
        channels = 1:size(unique_channels,2);
    end 
    
    n_channels = size(channels,2);
    
    [ begin_end, max_frequency] = extract_signal_set_data();

    sampling_ratio = SAMPLE_RESOLUTION / max_frequency;
    signal_duration =  begin_end(2) - begin_end(1);
    samples_in_signal = 1+ceil(signal_duration / sampling_ratio);    

    %If callback function is given, return initialize array large enough to
    %hold samples for a single channel and call the function with the
    %extracted header. Otherwise initialize an array for all channels.
    if is_function(callback)
        records = zeros(1, samples_in_signal);
        callback([],header,[]);
    else
        records = zeros(n_channels, samples_in_signal);       
    end

    for i = 1:size(unique_channels,2) 
                
        if ~ismember(unique_channels(i), channels)
            continue;
        end
        
        channel_part_files = file_map(unique_channels(i));     
   
        for c=1:size(channel_part_files,2)
                      
            full_file = char(strcat(ncs_path ,filesep, channel_part_files(c)));
            %Extract selected range from the current NCS file
            if isunix 
               [timestamps, NumberOfValidSamples, samples, ncs_header]   = Nlx2MatCSC_v3(full_file, [1,0,0,1,1],1,extraction_mode,extraction_vector);
            else
               [timestamps,  NumberOfValidSamples, samples, ncs_header]   = Nlx2MatCSC(full_file, [1,0,0,1,1],1,extraction_mode,extraction_vector);       
            end
            %Fill the record array with samples
            sample_index = 1;
            while sample_index < size(samples,2)
                %Ignore all records with incomplete set of samples
                if NumberOfValidSamples(sample_index) == SAMPLES_PER_RECORD
                    start = calculate_sample_index(timestamps(sample_index) - begin_end(1));
                    stop = start + SAMPLES_PER_RECORD-1;
                    
                    if is_function(callback)
                        records(1,start:stop) = reshape(samples(:,sample_index),1,SAMPLES_PER_RECORD)*header.ADBitVolts;   
                    else
                        records(i,start:stop) = reshape(samples(:,sample_index),1,SAMPLES_PER_RECORD)*header.ADBitVolts;                
                    end
                end
                
                sample_index = sample_index +1;
                
            end
            
            clear timestamps  NumberOfValidSamples  samples ncs_header;

        end
        %Call the callback with the samples extracted from the current
        %channel and initialize the record array with zeros.
        if is_function(callback)
            callback(records,header,i);
            records = zeros(1, samples_in_signal);
        end
        
    end   
    
    %Calculates the position of a sample in the record array
    function index = calculate_sample_index(time_point)
        index = 1+round( time_point / (SAMPLE_RESOLUTION / max_frequency));
    end
  
    function [ begin_end, max_frequency] = extract_signal_set_data()
        channel_part_files = file_map(channels(1));     
  
        begin_end = [ Inf 0];
        max_frequency = 0;
             
        for i = 1:size(unique_channels,2) 

            if ~ismember(unique_channels(i), channels)
                continue;
            end

            channel_part_files = file_map(unique_channels(i));     
        
            for c=1:size(channel_part_files,2)

                full_file = char(strcat(ncs_path ,filesep, channel_part_files(c)));
                %Extract selected range from the current NCS file
                if isunix 
                    [ timestamps,  ncs_header]  = Nlx2MatCSC_v3(full_file, [1,0,0,0,0],1,extraction_mode,extraction_vector);
                else
                    [ timestamps, ncs_header]  = Nlx2MatCSC(full_file, [1,0,0,0,0],1,extraction_mode,extraction_vector);       
                end

                if begin_end(1) > timestamps(1)
                   begin_end(1) = timestamps(1); 
                end
                if begin_end(2) < timestamps(end)
                   begin_end(2) = timestamps(end); 
                end            

                ADBitVolts = regexp(ncs_header(17,:),'(?<=(-ADBitVolts\s))[\s.0-9]*$','match');    
                ADBitVolts = str2double(ADBitVolts{1});
                header.ADBitVolts = ADBitVolts;
                
                channel_index = find( channels == unique_channels(i));

                header.time_created = regexp(ncs_header(8,:),'(?<=(-TimeCreated\s))[\s0-9:\\/]*$','match');
                header.time_closed = regexp(ncs_header(9,:),'(?<=(-TimeClosed\s))[\s0-9:\\/]*$','match');    
                frequency = regexp(ncs_header(15,:),'(?<=(-SamplingFrequency\s))[\s0-9]*$','match');
                frequency = str2double(frequency{1});
                header.frequency(1,channel_index) = frequency;   

                label = regexp(ncs_header(18,:),'(?<=(-AcqEntName\s))[\s0-9:A-Za-z\\/\-]*$','match');
                header.label{channel_index} = char(label{1});   
                header.units{channel_index} = 'V';           

                if max_frequency < frequency
                   max_frequency = frequency;
                end
                clear timestamps
            end
        end
    end
    
    if is_function(callback)
       records = [];
    end
end

