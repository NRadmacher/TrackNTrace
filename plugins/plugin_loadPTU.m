function [plugin] = plugin_loadPTU()

%    -------------- Definition of plugin --------------

% Name of the component these options are for
name = 'Import PTU';

% Type of plugin.
% 1: Candidate detection
% 2: Spot refinement/fitting
% 3: Tracking
% 4: Postprocessing
% 5: Import
% 6: Export
type = 5;

% The functions this plugin implements
mainFunc =  @read_PTU;

% Description of output parameters
outParamDescription = {'FLIM movie','head'}; % set in init function

% Create the plugin
plugin = TNTplugin(name, type, mainFunc, outParamDescription);

% Function to execute before frame-by-frame processing starts
% plugin.initFunc = @updateOutParamDescription;

% Description of plugin, supports sprintf format specifier like '\n' for a newline
plugin.info = struct(...
    'Description','Load a PTU file.',...
    'supportedFormats',{{'*.ptu','PTU'}},...
    'hasFLIM',true,...
    'hasTCSPC',true,...
    'getTCSPC',@getTCSPC);

% If hasTCPSC is true the plugin needs an accumulate function to retrieve
% the TCSPC of the localisations. This is given defined in plugin.postFunc.
% If plugin.initFunc is set the plugin supports caching via the output 
% 'index', which is passed instead of the inputfile to subsequent calls of 
% plugin.postFunc.
% The syntax has to be:
% PTU_accumulate(inputfile, outputs, framebinning, pixelbinning, timegate)
plugin.initFunc = @(opt)@PTU_accumulate;
plugin.postFunc = @accumulate_select;

plugin.add_param('alignBidirectional',...
    'bool',...
    false,...
    'Automatically determines and corrects the shift between lines for bidirectional scans.');
plugin.add_param('calibrationFile',...
    'filechooser',...
    {'','mat'},...
    'Select a file with ISM calibration data.');
plugin.add_param('fourierReweighting',...
    'bool',...
    false,...
    'Applies Fourier Reweighting to ISM shifted images to achieve resolution doubeling. If no calibrationFile is selected Fourier Reweighting is turned off');
plugin.add_param('ISM',...
    'bool',...
    true,...
    'Applies ISM shift vectors to photon positions. If no calibrationFile is selected Fourier Reweighting is turned off');
plugin.add_param('fastLT',...
    'list',...
    {'Std','Mean','Median','None','ChanRatio'},...
    'Chose how the fast lifetime is calculated. For ChanRatio, the relative intensity of the first channel is calculated.');
plugin.add_param('cacheMovie',...
    'int',...
    {2,0,10},...
    'Number of generated movie versions to save temporary. 0 disables cache. Oldest will be replaced if necessary.');

end


%   -------------- User functions --------------


function [movie,metadata] = read_PTU(pluginOptions,filename_movie, frame_range, frame_binning,fastLT,timegate)
%     profile on
    if endsWith(filename_movie,'ptu')
        if nargin<4 || isempty(frame_binning)
            frame_binning = 1;
        end
        if nargin<5 || isempty(fastLT) || contains(pluginOptions.fastLT,'none','IgnoreCase',true)
            fastLT = false;
        end
        if nargin<6
            timegate = [];
        end

        if isempty(pluginOptions.calibrationFile)
            % turn off fourierReweighting and ISM when no calibration file
            % is selected
            pluginOptions.fourierReweighting = false;
            pluginOptions.ISM = false;
            channelmap = [];
            shiftVector = [];
        else
            calib = open(pluginOptions.calibrationFile);
            n_pix = size(calib.shiftVector,2);
            channelmap = ones(n_pix,1);
            if pluginOptions.ISM
                shiftVector.sv          = calib.shiftVector;
                shiftVector.svPixelSize = calib.pixSize;
            else %confocal
                shiftVector = [];   
            end
        end
        if pluginOptions.fourierReweighting
            PSF = calib.psf;
            eps = 0.09;
        else
            PSF = [];
            eps = [];
        end

        % Check cache
        movieArgs = {filename_movie, frame_range, frame_binning,fastLT,timegate,pluginOptions.fastLT, shiftVector, PSF};
        outArgs = getCache(movieArgs{:});
        if ~isempty(outArgs)
            [movie,metadata] = outArgs{1:2};
            return;
        end

        if ~logical(fastLT)
            [head,movie] = PTU_accumulate(filename_movie,{'head','tag'},[frame_binning, frame_range],[],timegate,pluginOptions.alignBidirectional,channelmap,shiftVector);
            tau = [];
        else
            if ~isfield(pluginOptions,'fastLT')||isempty(pluginOptions.fastLT)
                pluginOptions.fastLT = 'Std';
            end
            switch pluginOptions.fastLT
                case 'Std'
                    [head,movie,tau] = PTU_accumulate(filename_movie,{'head','tag','tau'},[frame_binning, frame_range],[],timegate,pluginOptions.alignBidirectional,channelmap,shiftVector);
                case 'Mean'
                    [head,movie,tau] = PTU_accumulate(filename_movie,{'head','tag',{'tau',@mean}},[frame_binning, frame_range],[],timegate,pluginOptions.alignBidirectional,channelmap,shiftVector);
                case 'Median'
                    [head,movie,tau] = PTU_accumulate(filename_movie,{'head','tag',{'tau',@median}},[frame_binning, frame_range],[],timegate,pluginOptions.alignBidirectional,channelmap,shiftVector);
                case 'ChanRatio'
                    % create channelmap. We don't know yet how many
                    % channels are present. Assuming at most 64.
                    channelmap = [1;2*ones(63,1)];
                    [head,movie] = PTU_accumulate(filename_movie,{'head','tag'},[frame_binning, frame_range],[],timegate,pluginOptions.alignBidirectional,channelmap,shiftVector);
                    tau = movie(:,:,1,:); % tau is here the relative intensity in the first channel
                    movie = sum(movie,3); % The PTU channels are summed for now. In future this could be an option.
                    ind = tau>0;
                    tau(ind) = tau(ind)./movie(ind);
                    if head.TNTnChan>size(channelmap,1)
                        fprintf('Channel with number > 64 are ignored.');
                    end
                otherwise
                    error('Unknown fastLT parameter');
            end
            pluginOptions.uniqChan = head.TNTuniqChan; 
            tau = mean(permute(tau,[1 2 4 3]),4,"omitnan");
        end
        
        movie = sum(permute(movie,[1 2 4 3]),4); % The PTU channels are summed for now. In future this could be an option.
        if pluginOptions.fourierReweighting
            if calib.pixSize ~= head.TNTpixelSize
                %calculate psf for current pixel size
                PSF = calib.PSFfunc(calib.NA,calib.fd,calib.lamex,...
                    head.TNTpixelSize,calib.over);
            end
            fprintf('Fourier reweighting ... ');
            parfor i = 1:size(movie,3)
                movie(:,:,i) = ISM_frw(movie(:,:,i), PSF, eps);
            end
            fprintf('Done!\n');
        end

        if ~isempty(tau)
            movie = {movie,tau};
        else
            movie = {movie};
        end

        if nargout>1 || (isfield(pluginOptions,'cacheMovie') && pluginOptions.cacheMovie>0)
            metadata = struct(...
                'filename',filename_movie,...
                'pixelsize',head.ImgHdr_PixResol,...
                'pixelsize_unit',[char(181) 'm'],...% um
                'tau_unit','ns',...
                'bidirectional',logical(head.ImgHdr_BiDirect),...
                'ISM', ~isempty(pluginOptions.calibrationFile),...
                'FRW', pluginOptions.fourierReweighting,...
                'epsilon', eps,...
                'head',head...
                );
            if isfield(head,'ImgHdr_FrameFrequency')
                metadata.framerate = head.ImgHdr_FrameFrequency/frame_binning;
                metadata.framerate_unit = '1/s';
            end
        end
        % Save cache
        if isfield(pluginOptions,'cacheMovie')
            try
                setCache({movie,metadata},pluginOptions.cacheMovie, movieArgs{:});
            catch err
                warning('Could not cache the movie: Path might be read-only or contain special characters.');
                disp( getReport( err, 'extended', 'hyperlinks', 'on' ) );
            end
        end
    else
        error('No PTU file.');
    end
%     profile viewer
end

function [tcspcdata,resolution] = getTCSPC(inputfile,maxPhotons,pluginOptions)
    if endsWith(inputfile,'.mat')
        mfile = inputfile;
    else
        mfile = [inputfile(1:end-4) '_index.mat'];
    end
    
    if ~exist(mfile,'file')
        PTU_index(inputfile, pluginOptions.alignBidirectional);
        if ~exist(mfile,'file') % incase PTU_index fails
            error('Could not find PTU index file.');
        end
    end
    
    mf = matfile(mfile);
    head = mf.head;
    resolution = head.MeasDesc_Resolution; % in s
    
    tcspcdata = mf.im_tcspc(1:min(maxPhotons,size(mf,'im_tcspc',1)),1);
end

function argout = getCache(filename_movie, varargin)
    % varargin = {frame_range, frame_binning,fastLT,timegate,pluginOptions.fastLT, shiftVector, PSF}
    argout = {};
    try
        mf = openCache(filename_movie,false);
        if isempty(mf)
            return;
        end
        if any(strcmp(who(mf),'argout'))
            ncache = size(mf,'argout',1);
        else
            ncache = 0;
        end
        for icache = 1:ncache
            argin =  mf.argin(icache,1);
            argin = argin{1};
            fastLT = varargin{3};
            if ~fastLT
                % Extra FLIM data does not interfere
                varargin{5} = argin{6};
                varargin{3} = argin{4};
            end
            if isequaln(argin,{filename_movie,varargin{:}}) %#ok<CCAT> This is not identical with cat, since cat removes empty cells.
                % All arguments are equal.
                argout = mf.argout(icache,1);
                argout = argout{1};
                mf.time_read(icache,1) = now; % Update last read timestamp
                % Drop FLIM data if present and not requested
                if ~fastLT && numel(argout{1})>1
                    argout{1} = argout{1}(1);
                end
                break;
            end
        end
    catch err %#ok<NASGU>
        warning('Invalid cache file. Clearing cache.');
        openCache(filename_movie,false,true);  
    end
end

function setCache(argout,maxcache, filename_movie, varargin)
% varargin = {frame_range, frame_binning,fastLT,timegate,pluginOptions.fastLT, shiftVector, PSF}
    if maxcache == 0
        openCache(filename_movie,false,true);       
        return
    end
    mf = openCache(filename_movie,true);
    varlist = who(mf);
    if any(strcmp(varlist,'argout'))
        ncache = size(mf,'argout',1);
    else
        ncache = 0;
    end
    if ncache < maxcache 
        icache = ncache+1;
    else
        % find longest unread entry
        tread = mf.time_read;
        [~,icache] = min(tread);
    end
    mf.argin(icache,1) = {{filename_movie,varargin{:}}}; %#ok<CCAT> This is not identical with cat, since cat removes empty cells.
    mf.argout(icache,1) = {argout};
    mf.time_read(icache,1) = now;
    mf.time_created(icache,1) = now;
end

function mf = openCache(filename_movie,createFlag,clearFlag)
    cacheSuffix = '_FLIM.cache';
    [fpath,fname] = fileparts(filename_movie);
    cachename = fullfile(fpath,makeValidMatfileName([fname cacheSuffix]));% Matlab cannot write to matfiles with special characters in the name
    if exist(cachename,'file') && nargin > 2 && clearFlag
        delete(cachename);
    end
    if exist(cachename,'file') 
        % creates a matfile if it does not exist.
        mf = matfile(cachename,'Writable',true);        
    elseif (nargin>1 && createFlag)
        mf = matfile(cachename);
        mf.filename_movie = filename_movie;
    else
        mf = [];
    end
end

function accumFun = accumulate_select(opt)
    % opt importOptions
    if ~isempty(opt.calibrationFile)
        calib = open(opt.calibrationFile);
        shiftVector.sv          = calib.shiftVector;
        shiftVector.svPixelSize = calib.pixSize;
        n_pix = size(calib.shiftVector,2);
        channelmap = ones(n_pix,1);
        accumFun = @accumulate_shifted;
    else
        accumFun = @PTU_accumulate;
    end
    
    function varargout = accumulate_shifted(varargin)
        varargin{7} = channelmap;
        varargin{8} = shiftVector;
        [varargout{1:nargout}] = PTU_accumulate(varargin{:});
    end
end

