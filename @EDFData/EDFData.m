classdef EDFData < handle
  % EDFData  Class that represents data in an EDF file.
  % 
  % This class can be used to read EDF files. It can also be used to
  % convert EDF files in the lossless compressed MEF format when the user
  % also has access to MefWriter.jar. This can be downloaded from the
  % IEEG-Portal (www.ieeg.org).
  %
  % Author: J.B.Wagenaar
  % This code is loosely based on the EDF-Specification anf the EDFREAD
  % method made available on the Mathworks File Exchange.
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Copyright 2013 Trustees of the University of Pennsylvania
  % 
  % Licensed under the Apache License, Version 2.0 (the "License");
  % you may not use this file except in compliance with the License.
  % You may obtain a copy of the License at
  % 
  % http://www.apache.org/licenses/LICENSE-2.0
  % 
  % Unless required by applicable law or agreed to in writing, software
  % distributed under the License is distributed on an "AS IS" BASIS,
  % WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  % See the License for the specific language governing permissions and
  % limitations under the License.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  properties (SetAccess = private)
    ver = 0;
    patientID = '';
    recordID = '';
    startDate = '';
    startTime = '';
    bytes = 0;
    records = 0;
    duration = 0;
    ns = 0;
    label = {};
    transducer = '';
    units = {};
    prefilter = '';
    samples = [];
    sf = [];
    conversion = [];
    dcOffset = []
  end
  
  properties (SetAccess = private, Hidden)
    physicalMin = [];
    physicalMax = [];
    digitalMin = [];
    digitalMax = [];
    dataOffset
    dataPointer
    baseScale = [];
    simpleStruct = true;
    dpLabelStr = {};
    reserved = [];
    reserved2 = [];
  end
  
  methods
    function obj = EDFData(fileName)
      % EDFDATA  Represents an EDF file in Matlab.
      %
      %   OBJ = EDFDATA('filename') reads the specified EDF file and
      %   creates an object that can be used to read the data in the EDF
      %   file. 
      %
      %   The object contains the EDF header fields as properties and the
      %   datavalues can be accessed using the GETDATA and GETUNSCALEDDATA
      %   methods.
      %
      %   The CONVERT2MEF method converts the EDF file into a MEF file.
      %   This method can also be called with an array of EDFData files
      %   when multiple EDF files need to be concatenated.
      %
      
      [fid,msg] = fopen(fileName,'r');
      if fid == -1
        error(msg)
      end
      
      % HEADER
      obj.ver         = str2double(char(fread(fid,8)'));
      obj.patientID   = fread(fid,80,'*char')';
      obj.recordID    = fread(fid,80,'*char')';
      stDate          = fread(fid,8,'*char')';% (dd.mm.yy) or (dd-mm-yy)
      obj.startDate   = str2double(regexp(stDate,'\.|-','split'));
      stTime          = fread(fid,8,'*char')';% (hh.mm.ss) or (hh:mm:ss:)
      obj.startTime   = str2double(regexp(stTime,'\.|:','split'));
      obj.bytes       = str2double(fread(fid,8,'*char')');
      obj.reserved    = fread(fid,44);
      obj.records     = str2double(fread(fid,8,'*char')');
      obj.duration    = str2double(fread(fid,8,'*char')');
      
      % Number of signals
      obj.ns    = str2double(fread(fid,4,'*char')');
      for ii = 1:obj.ns
        obj.label{ii} = regexprep(fread(fid,16,'*char')','\W','');
      end
      % Transducer
      for ii = 1:obj.ns
        obj.transducer{ii} = fread(fid,80,'*char')';
      end
      % Physical dimension
      for ii = 1:obj.ns
        obj.units{ii} = fread(fid,8,'*char')';
      end
      % Physical minimum
      for ii = 1:obj.ns
        obj.physicalMin(ii) = str2double(fread(fid,8,'*char')');
      end
      % Physical maximum
      for ii = 1:obj.ns
        obj.physicalMax(ii) = str2double(fread(fid,8,'*char')');
      end
      % Digital minimum
      for ii = 1:obj.ns
        obj.digitalMin(ii) = str2double(fread(fid,8,'*char')');
      end
      % Digital maximum
      for ii = 1:obj.ns
        obj.digitalMax(ii) = str2double(fread(fid,8,'*char')');
      end
      for ii = 1:obj.ns
        obj.prefilter{ii} = fread(fid,80,'*char')';
      end
      for ii = 1:obj.ns
        obj.samples(ii) = str2double(fread(fid,8,'*char')');
      end
      for ii = 1:obj.ns
        obj.reserved    = fread(fid,32,'*char')';
      end  
      
      fclose(fid);
      
      % Set baseScale based on supplied units
      obj.baseScale = zeros(length(obj.units),1);
      unknownbaseScale = false;
      for ii = 1: length(obj.units)
        switch lower(deblank(obj.units{ii}))
          case {'v' 'volt' 'volts'}
            obj.baseScale(ii) = 1e6;
          case {'mv' 'millivolt' 'millivolts'}
            obj.baseScale(ii) = 1e3;
          case {'uv' 'microvolt' 'microvolts'}
            obj.baseScale(ii) = 1;
          case {'s' 'seconds' 'second'}
            obj.baseScale(ii) = 1;
          otherwise
            unknownbaseScale = true;
            obj.baseScale(ii) = 1; 
        end
      end
      if unknownbaseScale
        fprintf(['\nMethod could not determine the units of the data, '...
          'setting unknown conversion factors to 1.\n']);
      end
      
      % Find Conversion factor, DC-Offset, and sampling rates per channel.
      obj.conversion = (obj.physicalMax - obj.physicalMin)./(obj.digitalMax - obj.digitalMin);      
      obj.dcOffset = obj.physicalMax - obj.conversion .* obj.digitalMax;
      obj.sf = obj.samples./obj.duration;
      
      
      obj.dataOffset = 256 + obj.ns*256; % Fixed by EDF format.
      
      % If all sampling rates are the same, we can make a simple memmapfile
      % otherwise, we create a more complex memmapfile.
      if all(obj.samples == obj.samples(1))
        obj.simpleStruct = true;
        obj.dataPointer = memmapfile(fileName, ...
          'Offset',obj.dataOffset,...
          'Format', {'int16', [obj.samples(1) obj.ns obj.records],'x'});
      else
        obj.simpleStruct = false;
        names = regexp(sprintf('x%03i--',1:length(obj.label)),'--','split');
        names = names(1:end-1);
        obj.dpLabelStr = names;
        
        memFormat = cell(length(names),3);
        for iName = 1:length(names)
          memFormat(iName,:) = {'int16', [obj.samples(iName) 1], names{iName}};
        end
        
        obj.dataPointer = memmapfile(fileName, ...
          'Offset',obj.dataOffset, ...
          'Format', memFormat);  
          
        % Check size is the same as nrRecords
        if size(obj.dataPointer.Data,1) ~= obj.records
          warning('Nr records in header is not same as found');
          obj.records = size(obj.dataPointer.Data,1);
        end
      end
     
    end
    
    function out = getUnscaledData(obj, range, chIdx)
      % GETUNSCALEDDATA  Returns unscaled data
      %
      %   DATA = GETUNSCALEDDATA(OBJ, RANGE, CHIDX) returns the unscaled
      %   values of the EDF file, wher RANGE is a 1x2 vector of the start
      %   and end-index of the requested range and CHIDX is a vector of the
      %   channel indeces that are requested.
      
      assert(length(obj)==1,'OBJ should be a single EDFData object.');
      
      if obj.simpleStruct
        firstBlock = floor((range(1)-1)./obj.samples(1)) + 1; %ok
        firstOffset = mod(range(1)-1,obj.samples(1)); %ok

        lastBlock = ceil(range(2)./obj.samples(1)); %ok
        lastOffset = (lastBlock * obj.samples(1)) - range(2); %ok

        nblocks = lastBlock - firstBlock + 1;

        out = zeros(nblocks * obj.samples(1), length(chIdx));
        blockVec = firstBlock:lastBlock; 
        for iChan = 1: length(chIdx)
          aux = obj.dataPointer.Data.x(: ,chIdx(iChan), blockVec);
          out(:, iChan) = reshape(double(aux), numel(aux), 1);
        end
        
        % Chop Beginning end if necessary
        out = out((firstOffset+1):(length(out)-lastOffset),:);
        
      else
        
        out = zeros(range(2)-range(1)+1,length(chIdx));
        for iChan = 1: length(chIdx)

          firstBlock = floor((range(1)-1) ./obj.samples(chIdx(iChan))) + 1; %ok
          firstOffset = mod(range(1)-1,obj.samples(chIdx(iChan))); % ok
         
          lastBlock = ceil(range(2)./obj.samples(chIdx(iChan))); %ok
          lastOffset = (lastBlock * obj.samples(chIdx(iChan))) - range(2); %ok
          
          nblocks = lastBlock - firstBlock + 1;
          
          tempOut = zeros(nblocks * obj.samples(chIdx(iChan)), 1);
          curIdx = 1;
          for iBlock = firstBlock:lastBlock
            tempOut(curIdx:curIdx + obj.samples(chIdx(iChan))-1) = ...
              obj.dataPointer.Data(iBlock).(obj.dpLabelStr{chIdx(iChan)});
            curIdx = curIdx + obj.samples(chIdx(iChan));
          end
    
          % Chop Beginning end if necessary
          out(:,iChan) = tempOut((firstOffset+1):(length(tempOut)-lastOffset));
          
        end  
      end
    end
    
    function out = getData(obj, range, chIdx)
      % GETDATA  Returns scaled data
      %
      %   DATA = GETDATA(OBJ, RANGE, CHIDX) returns the data with the
      %   dc-offset and scaling factor applied. This data has the unit that
      %   is stored in the OBJ.Units property. RANGE is a 1x2 vector with
      %   the startIndex and endIndex of the requested subset of data.
      %   CHIDX is a 1xn vector of channel indeces.
      %
      %   For example:
      %     DATA = GETDATA(OBJ, [1 10000], 1:10)
      %
      %   See Also: EDFData.GETUNSCALEDDATA

      assert(length(obj)==1, 'OBJ should be a single EDFData object.');
      
      out = getUnscaledData(obj, range, chIdx);
      for iChan = 1: length(chIdx)
        out(:, iChan) = obj.dcOffset(chIdx(iChan)) + ...
          out(:, iChan).* obj.conversion(chIdx(iChan));
      end

    end
    
    function success = convert2Mef(obj, varargin)
      % CONVERT2MEF  Converts all channels to MEF files.
      %
      %   SUCCESS = CONVERT2MEF(OBJ) converts data stored in the EDF files
      %   represented by OBJ to MEF files. MEF stores data as one file per
      %   channel and uses lossless compression to minimize filesize. The
      %   timestamps are copied from the EDF files unless an optional
      %   'startTime' is provided.
      %
      %   OBJ can be a single EDFDATA object or an array of EDFDATA object
      %   in case the data in the various EDF-files should be concatenated.
      %
      %   The LABEL property is used to determine the filenames of the
      %   MEF-files. 
      %
      %   SUCCESS = CONVERT2MEF(OBJ, 'altNames', ALTNAMES, ...) converts
      %   the EDF data and stores the data per channel in user-defined
      %   filenames. ALTNAMES should be a structure that has the original
      %   channel labels as properties, and the alternative filenames as
      %   values for the properties.
      %
      %     i.e. 
      %     allLabels = {OBJ.label}
      %     altNames = {'Ch1' 'Ch2' 'Ch3'}
      %     altNamesStruct = cell2struct(altNames, allLabels, 2)
      %
      %   SUCCESS = CONVERT2MEF(OBJ, 'noDC', ...) omits the EDF DC-offset
      %   in the conversion,
      %
      %   SUCCESS = CONVERT2MEF(OBJ, 'startTime', STARTTIME) sets the
      %   timestamp of the first sample of the MEF file to STARTTIME, where
      %   STARTTIME is a string containing the time in the format:
      %   "hh:mm:ss". and uses a 24 hour clock. In this case, the date will
      %   be set to 1-1-2000.
      %
      %   SUCCESS = CONVERT2MEF(OBJ, 'startDate', STARTDATE) sets the
      %   timestamp of the date 
      %
      %   SUCCESS = CONVERT2MEF(OBJ, 'outputFolder', PATH) saves the MEF
      %   files to the folder specified in PATH. If the folder does not
      %   exist, it will be created.
      %
      %   SUCCESS = CONVERT2MEF(OBJ, 'hdr', HDRSTRUCT) Not fully
      %   implemented. Can be used to set header fields in the MEF files.
      
      
      % Default values
      STEPSIZE = 100000;  % Number of datapoints per MEF write.
      DEFAULTBLOCKSIZE = 5000; % Number of datapoints per MEF Block.
      useDC = true; % Convert DC Offset.
      hdrStruct = struct(); % Default will be empty header.
      MEFUpscale = 1; %Unused but can be used to multiply value in MEF file.
      useAltLabels = false;
      outputFolder = ''; %Current folder
      curIdx = 1;
      timeOffset = 0;
      if nargin > 1
        while curIdx <= length(varargin)
          assert(ischar(varargin{curIdx}), ...
            'Attribute name should be a string');
          
          switch varargin{curIdx}
            case 'noDC'
              useDC = false;
              curIdx = curIdx + 1;
            case 'hdr'
              assert(isstruct(varargin{curIdx+1}),...
                '''hdr'' parameter should be followed by struct.');
              hdrStruct = varargin{curIdx+1};
              curIdx = curIdx + 2;
            case 'outputFolder'
              assert(ischar(varargin{curIdx+1}),...
                'OutputFolder variable must be string.');
              outputFolder = varargin{curIdx+1};
              
              % Test exist folder, if not create
              if exist(outputFolder,'dir')~=7
                mkdir(outputFolder);
              end
              
              curIdx=curIdx+2;
            case 'altNames'
              assert(isstruct(varargin{curIdx+1}),...
                '''altNames'' parameter value should be structure.');
              altLabels = varargin{curIdx+1};
              allLabels = findAllLabels(obj);
              allPresent = all(cellfun(@(x) any(strcmp(x,allLabels)),fieldnames(altLabels)));
              assert(allPresent,'Not all Labels have an alternative name.');

              useAltLabels = true;
              curIdx = curIdx + 2;
            case 'startTime'
              assert(ischar(varargin{curIdx+1}),...
                '''startTime'' parameter must be of type string (hh:mm:ss)');
              try 
                % Get time offset in micro-seconds.
                timeOffset = str2double(regexp(varargin{curIdx+1},':','split'))*[3600 60 1]'*1e6;
              catch
                error('''StartTime'' parameter must be in format ''hh:mm:ss''');
              end
              curIdx = curIdx + 2;

            otherwise
             error('Incorrect attribute name');
          end
        end
      end

      % Find all available labels, and startTimes.
      for iObject = 1 : length(obj)
        tmpObject = obj(iObject);
        assert(~any(cellfun('isempty', tmpObject.label)),...
          'Some or all labels are empty, please set a label for each channel.');
        assert(length(unique(tmpObject.label)) == length(tmpObject.label),...
          'Labels within object are not unique, please reassign labels.');
      end
      
      allLabels = findAllLabels(obj);
      allStartTimes = zeros(length(obj),1);
      for iObject = 1 : length(obj)
        tmpObject = obj(iObject);
        allStartTimes(iObject) = 1e3*EDFData.datetime2utc(tmpObject.startDate(3),...
            tmpObject.startDate(2),tmpObject.startDate(1), tmpObject.startTime(1),...
            tmpObject.startTime(2),tmpObject.startTime(3));
      end      
      
      % sort concatObjects by StartTime.
      [allStartTimes, sortIx] = sort(allStartTimes);
      obj = obj(sortIx);
      
      % Set offsets for startTimes in case offsetTime is provided.
      firstObjOffset = allStartTimes(1) - 946684800000 + timeOffset;
      allStartTimes = allStartTimes - firstObjOffset;
      
      % Iterate over all unique labels
      for iLabel = 1: length(allLabels)
        display(sprintf('Converting Channel: %s',allLabels{iLabel}));   
        
        if useAltLabels
          curAltLabel = altLabels.(allLabels{iLabel});
          fileName = sprintf('%s.mef',curAltLabel);
          filePath = fullfile(outputFolder,fileName);
        else
          fileName = sprintf('%s.mef',allLabels{iLabel});
          filePath = fullfile(outputFolder,fileName);
        end
        
        curLabel = allLabels{iLabel};
        
        % Iterate over all files
        firstValidObject = true;
        for iObject = 1: length(obj)
          curObject = obj(iObject);
          curLabelIdx = find(strcmp(curLabel, curObject.label),1);
          
          % If label does not exist in File, continue to next file.
          if isempty(curLabelIdx)
            fprintf('-\n');
            continue;
          end
          
          % Find number of mef-writing steps for current file.
          nSteps = ceil((curObject.records * ...
            curObject.samples(curLabelIdx))./STEPSIZE);

          % If first file, creat mef-object
          if firstValidObject
            firstValidObject = false;
            curSF = curObject.sf(curLabelIdx);
            curBlockSize = ceil(DEFAULTBLOCKSIZE ./ curObject.sf(curLabelIdx)); 
            curGapThres = 1.25 * 1e6 * (1./curSF);
            
            mw = edu.mayo.msel.mefwriter.MefWriter(filePath, ...
              curBlockSize, curSF, curGapThres);

            % Find conversion factor to native Unit (uV in case of Volts).
            conv2uV = curObject.conversion(curLabelIdx);

            % Set MEF-Scale
            curOffset = curObject.dcOffset(curLabelIdx);
            curConv = curObject.conversion(curLabelIdx);
            curBase = curObject.baseScale(curLabelIdx);
            mw.setVoltageConversionFactor( (curConv * curBase) / MEFUpscale);
            
            % Determine ScaledOffset
            if useDC
              scaledOffset = curOffset / curConv;
            else
              scaledOffset = 0;
            end
            
            % Set other header info:
            hdrNames = fieldnames(hdrStruct);
            if ~isempty(hdrNames)
              for i = 1:length(hdrNames)
                switch hdrNames{i}
                  case 'subjectID'
                    mw.setSubjectID(hdrStruct(hdrNames{i}));
                  case 'institution'
                    mw.setInstitution(hdrStruct(hdrNames{i}));
                  otherwise
                    display(sprintf('Unknown header field: %s', hdrNames{i}));
                end
              end
            end
            
          elseif curSF ~= curObject.sf(curLabelIdx)
            fprintf(2,'Sampling Frequency mismatch for label: %s. Skipping file.',curLabel);
            continue;
          elseif curOffset ~= curObject.dcOffset(curLabelIdx);
            fprintf(2,'DC-Offset mismatch for label: %s. Skipping file.',curLabel);
            continue;
          elseif conv2uV ~= curObject.conversion(curLabelIdx);
            fprintf(2,'Conversion mismatch for label: %s. Skipping file.',curLabel);
            continue;
          end
          
          % Get StartTime in uUTC
          curTime = allStartTimes(iObject);
          timeStep = 1e6/curSF;

          curIdx = 1;
          progressIdx = ceil(linspace(1,nSteps,20)); % 5% increments in progressbar
          for iStep = 1:nSteps

            % Display progressbar
            if any(iStep==progressIdx)
              fprintf('.');
            end

            % Find current time vector
            timeVec = curTime : timeStep : (curTime + timeStep * STEPSIZE-1);

            % Get data for current step; if last step, get truncated data.
            if iStep < nSteps
              curData = scaledOffset + double(curObject.getUnscaledData(...
                [curIdx (curIdx + STEPSIZE - 1)], curLabelIdx));
            else
              totalLength = curObject.records * curObject.samples(curLabelIdx);
              curData =  scaledOffset +  double(curObject.getUnscaledData(...
                [curIdx totalLength], curLabelIdx));
               % Convert to baseUnit (uV)
              timeVec = timeVec(1 : length(curData));
            end

            curData = int32(MEFUpscale * curData);
            mw.writeData(curData, timeVec, length(curData));
            curTime = timeVec(end) + timeStep;
            curIdx= curIdx + STEPSIZE;
          end
          fprintf('\n');
          

        end
        mw.close();
      end
    success = true;
    end
    
    function allLabels = findAllLabels(obj)
      %FINDALLLABELS returns vector of all labels in supplied objects.
      %
      %   ALLLABELS = FINDALLLABELS(OBJ) returns a list of unique labels
      %   across the multiple suplied objects. OBJ is a single EDFData
      %   object or an array of EDFData objects.
      %
      %   If OBJ is a single EDFData object, this method returns the
      %   OBJ.LABEL property value as long as each label is unique.
      
      allLabels = {};
      for iObject = 1 : length(obj)
        allLabels = [allLabels obj(iObject).label]; %#ok<AGROW>
      end
        
      allLabels = unique(allLabels);
    end
    
    function obj = setLabel(obj, labelVec)
      % SETLABEL sets label names.
      %
      %   OBJ = SETLABEL(OBJ, LABELVEC) temporarily sets the labels in the
      %   EDFDATA object. This does not update the labels in the original
      %   EDF-File. This method can be used to provide labels for the
      %   EDF2MEF conversion method.
      
      assert(iscell(labelVec),'LabelVec must be cell-array.');
      assert(length(labelVec) == length(obj.label),...
        'Length of LabelVec must equal number of labels in object.');
      assert(all(cellfun(@(x) ischar(x),labelVec)),...
        'All labels must be strings');
      
      obj.label = labelVec;
      
      
    end
  end
  
  methods (Static)
    function out = datetime2utc(year, month, day, hour, minute, sec)
      %DATE2UTC  Converts date and duration to utc timestamps
      %
      % TIMESTAMP = DATE2UTC(YEAR,MONTH,DAY)
      %
      % [TIME1 TIME2] = DATE2UTC(YEAR, MONTH, DAY, SAMPLINGFREQ, NRSAMPLES)
      % provides two timestamps indicating the start and end-time of the
      % dataset in utc. 
      %
      % UTC Time is defined as the number of miliseconds from 1/1/1970. This
      % method does not take the time of the day into account, and therefor
      % always returns the start-time (TIME1/TIMESTAMP) as 12:00AM on that
      % particular day.
      %
      % By: Joost Wagenaar

      % Date
      if year < 1900
        year = year + 2000;
      end
      
      out = (datenum(year,month,day)- datenum(1970,1,1)) * 86400000;

      % Time
      out =out + (hour*3600 +minute*60 +sec)*1e3;

    end
    
    function out = find_edf_files(varargin)
      % FIND_EDF_FILES  Returns all EDF file-names in supplied path.
      %
      %   FILES = FIND_EDF_FILES returns EDF filenames in current path.
      %
      %   FILES = FIND_EDF_FILES(PATH) returns EDF files in supplied PATH.
      
      if nargin
        path = varargin{1};
      else
        path = pwd;
      end
      
      allFiles = dir(path);
      allFiles = allFiles(~[allFiles.isdir]);
      allFileNames = {allFiles.name};
      edfFiles = strcmp('.edf',cellfun(@(x) x(end-3:end),allFileNames,'UniformOutput',false));
      hiddenFiles = strcmp('.',cellfun(@(x) x(1),allFileNames,'UniformOutput',false));
      allFileNames = allFileNames(edfFiles & ~hiddenFiles);
      out = fullfile(path, allFileNames);

    end
        
  end
  
  
end