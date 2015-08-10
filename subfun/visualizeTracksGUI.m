function [h_main, selected_option] = visualizeTracksGUI(movie, trajectoryData, FPS, traj_lifetime, n_colors, use_bw, show_RunAgain)
% USAGE: visualizeTracksGUI(movie, trajectoryData)
% [ Full USAGE: visualizeTracksGUI(movie, trajectoryData, FPS, traj_lifetime, n_colors, use_bw) ]
%
% Visualizer for tracks computed by the tracker.
%
% Input:
%   movie: 3D matrix (rows,cols,frames) of the analyzed movie
%   trajectoryData: 2D matrix with columns particleID|frame|x|y|...
%                   This is the output of the tracker
% Inputs also adjustable by GUI:
%   FPS: frames per second to play movie with | default: 30
%   traj_lifetime: trajectories are kept for #traj_lifetime frames after
%                  the particle has vanished. | default: 0
%   n_colors: number of colors used to display trajectories | default: 20
%             Colors are generated using distinguishable_colors.m by
%             Timothy E. Holy (Matlab File Exchange).
%   use_bw: black/white image, otherwise colormap hot | default false
%   show_RunAgain: Display run again dialog after closing. Used by tracker
%                  demo
% Output:
%   h_main - Handle to the GUI figure
%   selected_option - Used for tracker demo mode. Displays a dialog wether
%                     or not to rerun the tracker with changed settings.
%
% Author: Simon Christoph Stein
% E-Mail: scstein@phys.uni-goettingen.de
% Date: 2015
%

if nargin<2
    fprintf(' Need input (movie,trajectoryData)!\n');
    return
end
num_argin = nargin;
parse_inputs();

% -- Preparing the GUI --
h_main = openfig('visualizeTracksGUI_Layout.fig');
set(h_main,'handleVisibility','on'); % Make figure visible to Matlab (might not be the case)
set(h_main,'CloseRequestFcn',@onAppClose); % For cleanup
set(h_main,'Toolbar','figure'); % Toolbar needed for zooming
set(h_main, 'DoubleBuffer', 'on') % Helps against flickering

h_all = guihandles(h_main);
axes(h_all.axes);

% Text on top
set(h_all.toptext,'String',sprintf('frame = 1/%i',size(movie,3)));

% Buttons
set(h_all.but_play,'Callback',@playCallback);
set(h_all.but_contrast,'Callback',@contrastCallback);

% Slider
set(h_all.slider,'Value',1, 'Min',1,'Max',size(movie,3),'SliderStep',[1/size(movie,3) 1/size(movie,3)],'Callback', @sliderCallback);
hLstn = handle.listener(h_all.slider,'ActionEvent',@updateSlider); % For continous update of the shown slider value
set(h_all.slider,'ButtonDownFcn',@onSliderClick);

% Edit fields
set(h_all.edit_FPS,'String',sprintf('%i',FPS), 'Callback', @fpsCallback);
set(h_all.edit_lifetime,'String',sprintf('%i',traj_lifetime), 'Callback', @lifetimeCallback);
set(h_all.edit_colors,'String',sprintf('%i',n_colors), 'Callback', @colorCallback);

% Checkbox
set(h_all.cb_bw, 'Value', use_bw, 'Callback',@bwCallback);

% Timer -> this controls playing the movie
warning('off','MATLAB:TIMER:RATEPRECISION');
h_all.timer = timer(...
    'ExecutionMode', 'fixedDelay', ...   % Run timer repeatedly
    'Period', 1/FPS, ...                % Initial period is 1 sec.
    'TimerFcn', @timed_update, ...
    'StartFcn', @onTimerStart, ...
    'StopFcn',  @onTimerStop); % Specify callback
warning('on','MATLAB:TIMER:RATEPRECISION');



% -- Preparation of plotting --
% Convert data into cell array which is better for plotting
% Each cell saves frame|x|y for one track
if isempty(trajectoryData)
    id_tracks = [];
    n_tracks = 0;
else
    id_tracks = unique(trajectoryData(:,1));
    n_tracks = numel(id_tracks);
end

cell_traj = cell(n_tracks ,1);
cnt = 1;
for iTrack = 1:n_tracks
    cell_traj{iTrack} = trajectoryData( trajectoryData(:,1)== id_tracks(cnt) , 2:4);
    cnt = cnt+1;
end

% Create the color pool
if use_bw
    bg = {'k','w'}; % background color
else
    bg = {'r','w'};
end

track_colors = repmat( distinguishable_colors(n_colors, bg), ceil(n_tracks/n_colors) ,1);
track_colors = track_colors(1:n_tracks,:);


% Plot first frame
% Set x,y,color limits
xl = [0.5,size(movie,2)+0.5];
yl = [0.5,size(movie,1)+0.5];
firstImg = movie(:,:,1);
zl = [min(firstImg(:)), max(firstImg(:))];

plotFrame(1);

xlim(xl);
ylim(yl);
caxis(zl);

% Variables for playback
timePerFrame = 1/FPS;
elapsed_time = 0;
frame = 1;

% Wait for figure to close
if(show_RunAgain)
    uiwait(h_main);
end


% --- Nested Functions (mostly callbacks) ---

% The main function of the application. This plays the movie if the
% timer is running
    function timed_update(timer, event)
        frame = frame+1;
        if(frame >= size(movie,3))
            frame = size(movie,3);
            updateTopText();
            
            stop(h_all.timer);
        end
        set(h_all.slider,'Value',frame);
        updateTopText()
        
        % Skip frame if computer is too slow
        if elapsed_time > timePerFrame
            elapsed_time = elapsed_time - timePerFrame;
            return;
        end
        tic_start = tic;
        updateFrameDisplay();
        elapsed_time = elapsed_time + toc(tic_start)- timePerFrame;
    end

    function onTimerStart(timer, event)
        set(h_all.but_play,'String','Pause');
        elapsed_time = 0;
    end

    function onTimerStop(timer, event)
        set(h_all.but_play,'String','Play');
    end

% Used to display the current frame as selected by the 'frame' variable
% Also this sets and saves the axis states (e.g. for zooming);
    function updateFrameDisplay()
        xl = xlim;
        yl = ylim;
        
        plotFrame(frame);
        
        % save the axis limits in case the user zoomed
        xlim(xl);
        ylim(yl);
        caxis(zl);
    end

% Plots the frame with the input index
    function plotFrame(iF)
        imagesc(movie(:,:,iF)); axis image; colormap gray;
        if use_bw
            colormap gray;
        else
            colormap hot;
        end
        %     title(sprintf('Frame %i/%i',iF,size(movie,3)));
        
        % Draw the tracks of currently visible particles
        hold on;
        for iTr = 1:n_tracks
            if iF < cell_traj{iTr}(1,1) || iF > cell_traj{iTr}(end,1) + traj_lifetime % skip all particles not visible (any longer)
                continue
            end
            mask_toPlot = cell_traj{iTr}(:,1)<=iF;
            plot(cell_traj{iTr}(mask_toPlot, 2), cell_traj{iTr}(mask_toPlot, 3), '.--','Color',track_colors(iTr,:));
        end
        hold off;
    end

% Switch play/pause by button
    function playCallback(hObj, eventdata)
        if frame == size(movie,3)
            frame = 1;
        end
        if strcmp(get(h_all.timer, 'Running'), 'off')
            start(h_all.timer);
        else
            stop(h_all.timer);
        end
    end

% Stop playing, adjust contrast, continue
    function contrastCallback(hObj, eventdata)
        isTimerOn = strcmp(get(h_all.timer, 'Running'), 'on');
        if isTimerOn
            stop(h_all.timer);
        end
        
        %This is done to prevent a warning dialog about display range outside data range
        axes(h_all.axes);
        currImg = movie(:,:,frame);
        zImg = [min(currImg(:)), max(currImg(:))];
        zl = [max(zImg(1),zl(1)), min(zImg(2),zl(2))];
        caxis(zl);
        
        him = imcontrast;
        uiwait(him);
        zl = caxis;
        
        if isTimerOn
            start(h_all.timer);
        end
    end

% Switch black-white and hot display mode
    function bwCallback(hObj, eventdata)
        isTimerOn = strcmp(get(h_all.timer, 'Running'), 'on');
        if isTimerOn
            stop(h_all.timer);
        end
        
        use_bw = ~use_bw;
        
        % create colors
        if use_bw % background color
            bg = {'k','w'};
        else
            bg = {'r','w'};
        end
        
        track_colors = repmat( distinguishable_colors(n_colors, bg), ceil(n_tracks/n_colors) ,1);
        track_colors = track_colors(1:n_tracks,:);
        
        if isTimerOn
            start(h_all.timer);
        else
            updateFrameDisplay()
        end
        
    end

% Update the movie FPS
    function fpsCallback(hObj, eventdata)
        isTimerOn = strcmp(get(h_all.timer, 'Running'), 'on');
        if isTimerOn
            stop(h_all.timer);
        end
        
        FPS = str2num(get(h_all.edit_FPS, 'String'));
        if isempty(FPS) || FPS<=0
            FPS = 30;
        end
        
        % Timer is limited to 1ms precision
        if FPS>1000
            FPS = 1000;
            warning('Max FPS is 1000 due to timer precision');
        end
        timePerFrame = 1/FPS;
        set(h_all.timer,'Period', timePerFrame);
        set(h_all.edit_FPS,'String',sprintf('%.1f',FPS));
        
        if isTimerOn
            start(h_all.timer);
        end
    end

% Update the lifetime of tracks
    function lifetimeCallback(hObj, eventdata)
        isTimerOn = strcmp(get(h_all.timer, 'Running'), 'on');
        if isTimerOn
            stop(h_all.timer);
        end
        
        traj_lifetime = round(str2num(get(h_all.edit_lifetime,'String')));
        if traj_lifetime<=0 || isempty(traj_lifetime)
            traj_lifetime = 0;
        end
        set(h_all.edit_lifetime,'String',sprintf('%i%',traj_lifetime));
        
        if isTimerOn
            start(h_all.timer);
        else
            updateFrameDisplay();
        end
    end

% Recompute the colors
    function colorCallback(hObj, eventdata)
        isTimerOn = strcmp(get(h_all.timer, 'Running'), 'on');
        if isTimerOn
            stop(h_all.timer);
        end
        
        n_colors = round(str2num(get(h_all.edit_colors,'String')));
        if n_colors<=0 || isempty(n_colors)
            n_colors = 1;
        end
        set(h_all.edit_colors,'String',sprintf('%i%',n_colors));
        
        % create colors
        if use_bw % background color
            bg = {'k','w'};
        else
            bg = {'r','w'};
        end
        
        track_colors = repmat( distinguishable_colors(n_colors, bg), ceil(n_tracks/n_colors) ,1);
        track_colors = track_colors(1:n_tracks,:);
        
        if isTimerOn
            start(h_all.timer);
        else
            updateFrameDisplay();
        end
    end

% This is called after letting the slider go
    function sliderCallback(hObj, eventdata)
        %
    end

% This is called continously when dragging the slider
    function updateSlider(hObj,eventdata)
        % Round slider value, sync text
        frame = round(get(h_all.slider,'Value'));
        set(h_all.slider,'Value',round(frame));
        updateTopText();
        
        % Stop timer
        if strcmp(get(h_all.timer, 'Running'), 'on')
            stop(h_all.timer);
        end
        
        % Skip frame if computer is too slow
        if elapsed_time >timePerFrame
            elapsed_time = elapsed_time - timePerFrame;
            return;
        end
        
        tic_start = tic;
        updateFrameDisplay();
        elapsed_time = elapsed_time + toc(tic_start)- timePerFrame;
    end

% Sets top text according to the current frame
    function updateTopText()
        set(h_all.toptext,'String',[sprintf('frame =  %i/%i',frame,size(movie,3))])
    end

% Parse input variables
    function parse_inputs()
        % input parsing
        if num_argin <3 || isempty(FPS)
            FPS = 30;
        end
        
        if num_argin <4 || isempty(traj_lifetime)
            traj_lifetime = 0;
        else
            traj_lifetime = round(traj_lifetime);
        end
        
        if num_argin < 5 || isempty(n_colors) || n_colors<=1
            n_colors = 20;
        else
            n_colors = round(n_colors);
        end
        
        if num_argin <6 || isempty(use_bw)
            use_bw = false;
        end
        
        if num_argin < 7 || isempty(show_RunAgain)
            show_RunAgain = false;
            selected_option = 'No';
        end
        
    end

% Cleanup function. This is neccessary to delete the timer!
    function onAppClose(hObj, event)
        if strcmp(get(h_all.timer, 'Running'), 'on')
            stop(h_all.timer);
        end
        delete(h_all.timer);
        
        % Dialog used in DEMO mode for getting return values.
        if show_RunAgain
            d = dialog('Position',[300 300 220 100],'Name','Run again?','WindowStyle','normal');
            
            txt = uicontrol('Parent',d,...
                'Style','text',...
                'Position',[5 40 210 40],...
                'String',sprintf('Run again to adjust settings?'));
            
            btn_yes = uicontrol('Parent',d,...
                'Position',[30 10 70 25],...
                'String','Yes',...
                'Callback',@buttonPress);
            
            btn_no = uicontrol('Parent',d,...
                'Position',[120 10 70 25],...
                'String','No',...
                'Callback',@buttonPress);
            selected_option = 'No';
            uiwait(d);
        end
        delete(h_main);
        
        function buttonPress(hObj, event)
            selected_option = get(hObj,'String');
            delete(d);
        end
    end

end
