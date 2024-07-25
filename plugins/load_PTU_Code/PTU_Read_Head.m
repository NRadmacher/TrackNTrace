function [head] = PTU_Read_Head(name)
    % Read PicoQuant Unified TTTR Files
    
    head = [];
    fid = fopen(name);
    if fid<1
        fprintf(1,'\n\n      Could not open <%s>. Aborted.\n', name);
    else
        
        tyEmpty8      = hex2dec('FFFF0008');
        tyBool8       = hex2dec('00000008');
        tyInt8        = hex2dec('10000008');
        tyBitSet64    = hex2dec('11000008');
        tyColor8      = hex2dec('12000008');
        tyFloat8      = hex2dec('20000008');
        tyTDateTime   = hex2dec('21000008');
        tyFloat8Array = hex2dec('2001FFFF');
        tyAnsiString  = hex2dec('4001FFFF');
        tyWideString  = hex2dec('4002FFFF');
        tyBinaryBlob  = hex2dec('FFFFFFFF');
            
        head.Magic =  deblank(char(fread(fid, 8, 'char')'));
        if not(strcmp(head.Magic,'PQTTTR'))
            error('Magic invalid, this is not an PTU file.');
        end
        head.Version = deblank(char(fread(fid, 8, 'char')'));
        
        
        TagIdent = deblank(char(fread(fid, 32, 'char')'));  % TagHead.Ident
        TagIdx   = fread(fid, 1, 'int32');                  % TagHead.Idx
        TagTyp   = fread(fid, 1, 'uint32');                 % TagHead.Typ
        
        if TagIdent(1)==36
            TagIdent = ['Remote_' TagIdent(2:end)];
        end
        
        while ~strcmp(TagIdent, 'Header_End')
                    
            if TagIdx > -1
                EvalName = ['head.' TagIdent '(' int2str(TagIdx + 1) ')'];
            else
                EvalName = ['head.' TagIdent];
            end
            
            % check Typ of Header
            
            switch TagTyp
                case tyEmpty8
                    fread(fid, 1, 'int64');
                    eval([EvalName '= [];']);
                case tyBool8
                    TagInt = fread(fid, 1, 'int64');
                    if TagInt==0
                        eval([EvalName '=false;']);
                    else
                        eval([EvalName '=true;']);
                    end
                case tyInt8
                    TagInt = fread(fid, 1, 'int64');
                    eval([EvalName '=TagInt;']);
                case tyBitSet64
                    TagInt = fread(fid, 1, 'int64');
                    eval([EvalName '=TagInt;']);
                case tyColor8
                    TagInt = fread(fid, 1, 'int64');
                    eval([EvalName '=TagInt;']);
                case tyFloat8
                    TagFloat = fread(fid, 1, 'double');
                    eval([EvalName '=TagFloat;']);
                case tyFloat8Array
                    TagInt = floor(fread(fid, 1, 'int64')/8);
                    TagArray = fread(fid, TagInt, 'double');
                    eval([EvalName '=TagArray;']);
                case tyTDateTime
                    TagFloat = fread(fid, 1, 'double');
                    eval([EvalName '=datestr(693960+TagFloat);']); 
                case tyAnsiString
                    TagInt = fread(fid, 1, 'int64');
                    TagString = deblank(char(fread(fid, TagInt, 'char'))');
                    tmp = char(regexp(EvalName,'\([0-9]+\)','match'));
                    if ~isempty(tmp)
                        EvalName=strrep(EvalName,tmp,['{' tmp(2:end-1) '}']);
                    end
                    eval([EvalName '=TagString;']);
                case tyWideString
                    % Matlab does not support Widestrings at all, just read and
                    % remove the 0's (up to current (2012))
                    TagInt = fread(fid, 1, 'int64');
                    TagString = fread(fid, TagInt, '*char');
                    TagString = (TagString(TagString ~= 0))';
                    fprintf(1, '%s\n', TagString);
                    if TagIdx > -1
                        EvalName = [TagIdent '(' int2str(TagIdx + 1) ',:)'];
                    end
                    %ignore due to luminosa header. usually not important
                    %eval([EvalName '=TagString;']);
                case tyBinaryBlob
                    TagInt = fread(fid, 1, 'int64');
                    TagBytes = fread(fid, TagInt, 'uint8');
                otherwise
                    
            end
            TagIdent = deblank(char(fread(fid, 32, 'char')'));  % TagHead.Ident
            TagIdx   = fread(fid, 1, 'int32');                  % TagHead.Idx
            TagTyp   = fread(fid, 1, 'uint32');                 % TagHead.Typ
        end
        
        head.length = ftell(fid)+8;
        fclose(fid);
    end
    %% Argument parser for consistant TNT naming convention
    
    %possible field names acros platforms
    dwellTimeFields     = {'ImgHdr_TimePerPixel'};
    pixelSizeFields     = {'ImgHdr_PixResol'};
    isBiDirFields       = {'ImgHdr_BiDirect'};
    frameRateFields     = {''};
    tcspsResFields      = {'MeasDesc_Resolution'};
    dimensionFields     = {'ImgHdr_Dimensions'};
    identFields         = {'ImgHdr_Ident'};
    pixXFields          = {'ImgHdr_PixX'};
    pixYFields          = {'ImgHdr_PixY'};
    syncRateFields      = {'TTResult_SyncRate'};
    numberOfChFields    = {'HW_InpChannels'};
    numberOfRecFields   = {'TTResult_NumberOfRecords'};

    
    %dwell time per pixel [s]
    for i = 1:numel(dwellTimeFields)
        if isfield(head,dwellTimeFields{i})
            head.TNTdwellTime = head.(dwellTimeFields{i});
        end
    end

    %pixel size [m] 
    for i = 1:numel(pixelSizeFields)
        if isfield(head,pixelSizeFields{i})
            head.TNTpixelSize = head.(pixelSizeFields{i});
        end
    end

    %is scanner running in both directions [logic]
    for i = 1:numel(isBiDirFields)
        if isfield(head,isBiDirFields{i})
            head.TNTisBiDirectional = head.(isBiDirFields{i});
        end
    end

    %TCSPC bin size [s] 
    for i = 1:numel(tcspsResFields)
        if isfield(head,tcspsResFields{i})
            head.TNTtcspsBinSize = head.(tcspsResFields{i});
        end
    end

    %dimension of the aquision [#] 
    for i = 1:numel(dimensionFields)
        if isfield(head,dimensionFields{i})
            head.TNTmeasDim = head.(dimensionFields{i});
        end
    end

    %aquision identifier[#] 
    for i = 1:numel(identFields)
        if isfield(head,identFields{i})
            head.TNTident = head.(identFields{i});
        end
    end

    %number of pixles in x direction [#] 
    for i = 1:numel(pixXFields)
        if isfield(head,pixXFields{i})
            head.TNTpixX = head.(pixXFields{i});
        end
    end

    %number of pixles in y direction [#] 
    for i = 1:numel(pixYFields)
        if isfield(head,pixYFields{i})
            head.TNTpixY = head.(pixYFields{i});
        end
    end

    %Laser sync rate [Hz] 
    for i = 1:numel(syncRateFields)
        if isfield(head,syncRateFields{i})
            head.TNTsyncRate = head.(syncRateFields{i});
        end
    end

    %number of channels or detectors [#] 
    for i = 1:numel(numberOfChFields)
        if isfield(head,numberOfChFields{i})
            %this filed contains the data channels
            if isfield(head,'TTResult_InputRate')%TODO are there other names for this channel
                %only use chanels that record photons
                if sum((head.TTResult_InputRate ~= 0)) < head.(numberOfChFields{i})
                    ind = find(head.TTResult_InputRate);
                    if ~isempty(ind)
                        %We asume that hydra harp is not used for ISM
                        if head.(numberOfChFields{i}) <= 8
                            head.TNTnChan       = numel(ind);
                            head.TNTchanMap     = ones(numel(head.TTResult_InputRate),1);
                        %in all other cases it is assumed no channels are
                        %skipped
                        else
                            head.TNTnChan       = ind(end);
                            head.TNTchanMap     = ones(head.TNTnChan,1);
                        end
                    else %if all Imputrates are 0. Currently fix for luminosa data
                        [~,~,ic] = unique(head.HWInpChan_TrgLevel);
                        ind = (ic.' == 1)&(head.HWInpChan_Offset ~= 0);
                        head.TNTnChan       = sum(ind);
                        head.TNTchanMap     = ones(head.TNTnChan,1);
                    end
                end
            else
                head.TNTnChan       = head.(numberOfChFields{i});
                head.TNTchanMap     = ones(head.TNTnChan,1);
            end
        end
    end

    %number of photons recorded [#] 
    for i = 1:numel(numberOfRecFields)
        if isfield(head,numberOfRecFields{i})
            head.TNTnPhoton = head.(numberOfRecFields{i});
        end
    end

end
