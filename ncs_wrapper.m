function [ records header ] = ncs_wrapper( ncs_path, sample_range, channels  )
%NCS_WRAPPER Wrapper for Neuralynx NCS format
%   Encapsulates the reading of Neuralunx NCS files. The name of the
%   channel recorded in the file is assumed to be presented in the last
%   numeric value of the file name, preceded by -, e.g. Cage1-1.ncs to
%   indicate channel 1 from cage 1.
%   Parameters:
%       ncs_path: Path to folder containing the NCS files.
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
%  	Neuralynx MATLAB Import/Export MEX Files for Windows http://neuralynx.com/software/NeuralynxMatlabImportExport_v6.0.0.zip
%	Neuralynx MATLAM Import/Export MEX Files for Linux http://neuralynx.com/software/Nlx2Mat_relDec15.tar.gz

    SAMPLES_PER_RECORD = 512;
    
    extraction_mode = 1;
    extraction_vector = 1;
    if nargin > 1 && ~isempty(sample_range)
        extraction_mode = 2;
        extraction_vector = [1 max(floor(sample_range(2)/SAMPLES_PER_RECORD),1)];       
    end
    
    records = [];   
    header = struct;    
    files = dir(fullfile(ncs_path, '*.ncs'));
     
    if nargin < 3
        channels = 1:size(files,1);
    end
    %Each file in the folder is assumed to present an unique channel
    for c = 1:size(files,1)

        full_file = strcat(ncs_path ,filesep, files(c).name);
        %Extract selected range from the current NCS file
        if isunix 
            [ samples, ncs_header]  = Nlx2MatCSC_v3(full_file, [0,0,0,0,1],1,extraction_mode,extraction_vector);
        else
            [ samples, ncs_header]  = Nlx2MatCSC(full_file, [0,0,0,0,1],1,extraction_mode,extraction_vector);       
        end
        n_samples = size(samples,2)*SAMPLES_PER_RECORD;
       
        %Channel number is assumed to be presented as the last numeric
        %value of the file name.
        channel_number = str2double(regexp(files(c).name, '(?<=[-])[\d]*(?=(.ncs))', 'match'));
        
        %Skip unselected channels        
        if ~ismember(channel_number,channels)
           continue 
        end            
        
        channel_index = find( channels == channel_number);

        if isempty(records)
           records = zeros(size(channels,2), n_samples ); 
           header.time_created = regexprep(ncs_header(8,:),'-TimeCreated ','');
           header.time_closed = regexprep(ncs_header(9,:),'-TimeClosed  ','');       
           header.sampling_frequency = regexprep(ncs_header(15,:),'-SamplingFrequency  ','');           
        end
        records(channel_index,:) = reshape(samples,1,n_samples);

    end
end


