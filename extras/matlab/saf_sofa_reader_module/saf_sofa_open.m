function s_sofa = saf_sofa_open(sofa_filepath)
% A "no-nonsense" function that uses C-style netcdf4 calls to extract SOFA
% variables, attributes and global attributes, from a .sofa file.
%
% sofa_filepath should include the .sofa file extention
%
% This function conforms to the SOFA 1.0 FIR convention defined in:
%     https://www.sofaconventions.org/mediawiki/index.php/GeneralFIR
%
% Author: Leo McCormack
% Date: 04.10.2020
% License: ISC

%% Info on sofa file
ncdisp(sofa_filepath)

%% Open sofa file and determine dimension IDs and lengths
ncid = netcdf.open(sofa_filepath, 'NOWRITE'); % read-only access
[numdims,numvars,numglobalatts,unlimdimid] = netcdf.inq(ncid); %#ok
for i=1:numdims 
    [dimname(i,1), dimlen(i,1)] = netcdf.inqDim(ncid,i-1); %#ok
    dimid(i,1) = netcdf.inqDimID(ncid,dimname(i,1));       %#ok
end  

%% Default Vars
s_sofa.Data_IR = [];
s_sofa.Data_SamplingRate = 0;
s_sofa.Data_Delay = 0;
s_sofa.SourcePosition = [];
s_sofa.ListenerPosition = []; 
s_sofa.ListenerUp = [];
s_sofa.ListenerView = [];
s_sofa.EmitterPosition = [];
s_sofa.ReceiverPosition = [];
s_sofa.MeasurementSourceAudioChannel = [];
s_sofa.MeasurementAudioLatency = [];

%% Default Attributes
s_sofa.ListenerPosition_Type = 'unknown';
s_sofa.ListenerPosition_Units = 'unknown';
s_sofa.ReceiverPosition_Type = 'unknown';
s_sofa.ReceiverPosition_Units = 'unknown';
s_sofa.SourcePosition_Type = 'unknown';
s_sofa.SourcePosition_Units = 'unknown';
s_sofa.EmitterPosition_Type = 'unknown';
s_sofa.EmitterPosition_Units = 'unknown';
s_sofa.Data_SamplingRate_Units = 'unknown';

%% Loop over vars
for i=1:numvars 
    [varnames{i},xtype,dimids,natts] = netcdf.inqVar(ncid,i-1); %#ok
    varname = varnames{i};
    switch varname 
        case 'Data.IR' 
            s_sofa.Data_IR = netcdf.getVar(ncid, netcdf.inqVarID(ncid,varname));  
        case 'Data.SamplingRate' 
            s_sofa.Data_SamplingRate = netcdf.getVar(ncid, netcdf.inqVarID(ncid,varname)); 
            for j=1:natts
                attname = netcdf.inqAttName(ncid,i-1,j-1);
                switch attname 
                    case 'Units'
                        s_sofa.Data_SamplingRate_Units = netcdf.getAtt(ncid,i-1,attname);
                end 
            end
        case 'Data.Delay'
            s_sofa.Data_Delay = netcdf.getVar(ncid, netcdf.inqVarID(ncid,varname));  
        case 'SourcePosition' 
            s_sofa.SourcePosition = netcdf.getVar(ncid, netcdf.inqVarID(ncid,varname)); 
            for j=1:natts
                attname = netcdf.inqAttName(ncid,i-1,j-1);
                switch attname
                    case 'Type'
                        s_sofa.SourcePosition_Type = netcdf.getAtt(ncid,i-1,attname);
                    case 'Units'
                        s_sofa.SourcePosition_Units = netcdf.getAtt(ncid,i-1,attname);
                end 
            end
        case 'ReceiverPosition'
            s_sofa.ReceiverPosition = netcdf.getVar(ncid, netcdf.inqVarID(ncid,varname)); 
            for j=1:natts
                attname = netcdf.inqAttName(ncid,i-1,j-1);
                switch attname
                    case 'Type'
                        s_sofa.ReceiverPosition_Type = netcdf.getAtt(ncid,i-1,attname);
                    case 'Units'
                        s_sofa.ReceiverPosition_Units = netcdf.getAtt(ncid,i-1,attname);
                end 
            end
        case 'ListenerPosition' 
            s_sofa.ListenerPosition = netcdf.getVar(ncid, netcdf.inqVarID(ncid,varname)); 
            for j=1:natts
                attname = netcdf.inqAttName(ncid,i-1,j-1);
                switch attname
                    case 'Type'
                        s_sofa.ListenerPosition_Type = netcdf.getAtt(ncid,i-1,attname);
                    case 'Units'
                        s_sofa.ListenerPosition_Units = netcdf.getAtt(ncid,i-1,attname);
                end 
            end
        case 'ListenerUp' 
            s_sofa.ListenerUp = netcdf.getVar(ncid, netcdf.inqVarID(ncid,varname)); 
        case 'ListenerView'
            s_sofa.ListenerView = netcdf.getVar(ncid, netcdf.inqVarID(ncid,varname)); 
        case 'EmitterPosition' 
            s_sofa.EmitterPosition = netcdf.getVar(ncid, netcdf.inqVarID(ncid,varname));  
            for j=1:natts
                attname = netcdf.inqAttName(ncid,i-1,j-1);
                switch attname
                    case 'Type'
                        s_sofa.EmitterPosition_Type = netcdf.getAtt(ncid,i-1,attname);
                    case 'Units'
                        s_sofa.EmitterPosition_Units = netcdf.getAtt(ncid,i-1,attname);
                end 
            end
        case 'MeasurementSourceAudioChannel'
            s_sofa.MeasurementSourceAudioChannel = netcdf.getVar(ncid, netcdf.inqVarID(ncid,varname));  
        case 'MeasurementAudioLatency'
            s_sofa.MeasurementAudioLatency = netcdf.getVar(ncid, netcdf.inqVarID(ncid,varname));  
    end
end

%% Default Global Atts
s_sofa.DataType = 'unknown';
s_sofa.Conventions = 'unknown';
s_sofa.Version = 'unknown';
s_sofa.SOFAConventions = 'unknown';
s_sofa.SOFAConventionsVersion = 'unknown';
s_sofa.APIName = 'unknown';
s_sofa.APIVersion = 'unknown';
s_sofa.ApplicationName = 'unknown';
s_sofa.ApplicationVersion = 'unknown';
s_sofa.Organization = 'unknown';
s_sofa.License = 'unknown';
s_sofa.DataType = 'unknown';
s_sofa.RoomType = 'unknown';
s_sofa.Origin = 'unknown';
s_sofa.DateCreated = 'unknown';
s_sofa.DateModified = 'unknown';
s_sofa.Title = 'unknown';
s_sofa.DatabaseName = 'unknown';
s_sofa.ListenerShortName = 'unknown';
s_sofa.AuthorContact = 'unknown';
s_sofa.History = 'unknown';
s_sofa.References = 'unknown';
s_sofa.Comment = 'unknown';
s_sofa.ReceiverDescription = 'unknown';
s_sofa.RoomDescription = 'unknown';
s_sofa.RoomLocation = 'unknown';
s_sofa.SourceDescription = 'unknown';
s_sofa.EmitterDescription = 'unknown';

%% Loop over Global attributes
for i=1:numglobalatts  
    attnames{i} = netcdf.inqAttName(ncid,netcdf.getConstant('NC_GLOBAL'),i-1);  
    attname = attnames{i};
    switch attname
        case 'DataType'
            s_sofa.DataType = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),attname);
        case 'Conventions'
            s_sofa.Conventions = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),attname);
        case 'Version' 
            s_sofa.Version = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),attname);
        case 'SOFAConventions'
            s_sofa.SOFAConventions = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),attname);
        case 'SOFAConventionsVersion'
            s_sofa.SOFAConventionsVersion = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),attname);
        case 'APIName'
            s_sofa.APIName = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),attname);
        case 'APIVersion'
            s_sofa.APIVersion = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),attname); 
        case 'ApplicationName'
            s_sofa.ApplicationName = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),attname);
        case 'ApplicationVersion'
            s_sofa.ApplicationVersion = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),attname);
        case 'Organization'
            s_sofa.Organization = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),attname);
        case 'License'
            s_sofa.License = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),attname);
        case 'RoomType'
            s_sofa.RoomType = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),attname);
        case 'Origin'
            s_sofa.Origin = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),attname);
        case 'DateCreated'
            s_sofa.DateCreated = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),attname);
        case 'DateModified'
            s_sofa.DateModified = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),attname);
        case 'Title'
            s_sofa.Title = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),attname);
        case 'DatabaseName'
            s_sofa.DatabaseName = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),attname);
        case 'ListenerShortName'
            s_sofa.ListenerShortName = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),attname);
        case 'AuthorContact'
            s_sofa.AuthorContact = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),attname);
        case 'History' 
            s_sofa.History = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),attname);
        case 'References'
            s_sofa.References = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),attname);
        case 'Comment'
            s_sofa.Comment = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),attname);
        case 'ReceiverDescription'
            s_sofa.ReceiverDescription = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),attname);
        case 'RoomLocation'
            s_sofa.RoomLocation = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),attname);
        case 'SourceDescription'
            s_sofa.SourceDescription = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),attname);
        case 'EmitterDescription'
            s_sofa.EmitterDescription = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),attname);
    end
end

netcdf.close(ncid);

end
