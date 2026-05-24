function varargout = GUI(varargin)
% GUI MATLAB code for GUI.fig
%      GUI, by itself, creates a new GUI or raises the existing
%      singleton*.
%
%      H = GUI returns the handle to a new GUI or the handle to
%      the existing singleton*.
%
%      GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI.M with the given input arguments.
%
%      GUI('Property','Value',...) creates a new GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI

% Last Modified by GUIDE v2.5 07-Feb-2018 11:51:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @GUI_OpeningFcn, ...
    'gui_OutputFcn',  @GUI_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    %  try
    gui_mainfcn(gui_State, varargin{:});
    %  catch me
    %      display(me)
    %  end
end
% End initialization code - DO NOT EDIT


% --- Executes just before GUI is made visible.
function GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI (see VARARGIN)

% Set up the data structures
if(length(varargin)==1)
    
    results = varargin{1};
    
    % Initialise menu handles
    set(handles.xmenu,'String',results.paramDefCellInput(:,1));
    set(handles.ymenu,'String',results.paramDefCellInput(:,1));
    
    addfomtosub = false;
    if sum(results.allEvaluationValues(1) == results.allSubEvaluationValues(:,1)) == 0
        addfomtosub = true;
        handles.fomat = 1;
    else
        [~,fomloc] = max(results.allEvaluationValues(1) == results.allSubEvaluationValues(:,1));
        handles.fomat = fomloc;
    end
    
    
    if isfield(results.DEParams,'subFoMNames') == 1
        
        names = results.DEParams.subFoMNames;
        if addfomtosub
            names = [{'FoM'},names];
        end
        set(handles.OptMenu,'String',names);
        
        set(handles.normalisation,'String',[{'1'},names]);
    else
        nosub = size(results.allSubEvaluationValues,1);
        if addfomtosub
            nosub = nosub +1;
        end
        
        names = cell(1, nosub);
        for i = 1:nosub
            names(i) = {strcat('Variable ',num2str(i))};
        end
        set(handles.OptMenu,'String',names);
        
        set(handles.normalisation,'String',[{'1'},names]);
    end
    
    
    
    % Get data for plotting
    handles.xydata=varargin{1}.allTestedMembers;
    
    if addfomtosub
        handles.zdata=[results.allEvaluationValues;results.allSubEvaluationValues];
    else
        handles.zdata=results.allSubEvaluationValues;
    end
elseif (length(varargin)==3)
    set(handles.xmenu,'String',varargin{1});
    set(handles.ymenu,'String',varargin{1});
    
    % Get data for plotting
    handles.xydata=varargin{2};

    handles.zdata=varargin{3};
    
elseif (length(varargin)==4)
    set(handles.xmenu,'String',varargin{1});
    set(handles.ymenu,'String',varargin{1});
    
    % Get data for plotting
    handles.xydata=varargin{2};

    addfomtosub = false;
    if sum(varargin{3}(1) == varargin{3}(:,1)) == 0
        addfomtosub = true;
    end
    
    if addfomtosub
        handles.zdata=[varargin{3};varargin{4}];
    else
        handles.zdata=varargin{4};
    end
    
else
    error('Data not entered into GUI correctly.')
end

% Set initial plots
handles.xnum=1;

contents = get(handles.xmenu,'string');
handles.xname=contents{get(handles.xmenu,'Value')};

handles.ynum=2;

contents = get(handles.ymenu,'string');
handles.yname=contents{get(handles.ymenu,'Value')};

contents = get(handles.OptMenu,'string');
handles.zname=contents{get(handles.OptMenu,'Value')};
if(get(handles.normalisation,'value')~=1)
    contents = get(handles.normalisation,'string');
    handles.zname = strcat(handles.zname,'/',contents{get(handles.normalisation,'Value')});
end

handles.datanum=1;
handles.normnum=1;
handles.zoomfrac=0;
handles.trimnegative = 0;
handles.trimpositive = 0;
handles.plotoutliers = 1;
handles.varno = 1;
handles.highlighttopFoM = 0;
handles.highlighttopresult = 0;

handles.nantoggle=0;


drawgraph(handles);

% Choose default command line output for GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in xmenu.
function xmenu_Callback(hObject, eventdata, handles)
% hObject    handle to xmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns xmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from xmenu

contents = cellstr(get(hObject,'String'));
handles.xnum = get(hObject,'Value');
handles.xname = contents{get(hObject,'Value')};
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function xmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ymenu.
function ymenu_Callback(hObject, eventdata, handles)
% hObject    handle to ymenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ymenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ymenu

contents = cellstr(get(hObject,'String'));
handles.ynum = get(hObject,'Value');
handles.yname = contents{get(hObject,'Value')};
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function ymenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ymenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plotbutton.
function plotbutton_Callback(hObject, eventdata, handles)
% hObject    handle to plotbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

drawgraph(handles);

guidata(hObject, handles);


% --- Executes on selection change in OptMenu.
function OptMenu_Callback(hObject, eventdata, handles)
% hObject    handle to OptMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns OptMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from OptMenu

handles.datanum = get(hObject,'Value');
contents = cellstr(get(hObject,'String'));
handles.zname = contents{get(hObject,'Value')};

if(get(handles.normalisation,'value')~=1)
    contents = get(handles.normalisation,'string');
    handles.zname = strcat(handles.zname,'/',contents{get(handles.normalisation,'Value')});
end

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function OptMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OptMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in normalisation.
function normalisation_Callback(hObject, eventdata, handles)
% hObject    handle to normalisation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns normalisation contents as cell array
%        contents{get(hObject,'Value')} returns selected item from normalisation
handles.normnum = get(hObject,'Value');

contents = cellstr(get(handles.OptMenu,'String'));
handles.zname = contents{get(handles.OptMenu,'Value')};

if(get(handles.normalisation,'value')~=1)
    contents = get(handles.normalisation,'string');
    handles.zname = strcat(handles.zname,'/',contents{get(handles.normalisation,'Value')});
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function normalisation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to normalisation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

value = get(hObject,'Value');
min = get(hObject,'Min');
max = get(hObject,'Max');
handles.zoomfrac = (value-min)/(max-min);

guidata(hObject, handles);

drawgraph(handles);

guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function drawgraph(handles)

%Initial plots

data=handles.zdata(handles.datanum,:);
fom = handles.zdata(handles.fomat,:);


if (handles.normnum>1)
    data = data./handles.zdata(handles.normnum-1,:);
end

% Choose what data to remove

test = zeros(1,length(data));
if handles.trimpositive
    test = test + (data - nanmean(data)) > handles.varno*nanstd(data);
end

if handles.trimnegative
    test = ((test+ (nanmean(data) - data) > handles.varno*nanstd(data))/2);
end

I = test > 0;

outliers = excludedata(handles.xydata(handles.ynum,:),data,'indices',I);

if strcmp(version,'8.1.0.604 (R2013a)')
    mid = nanmean(data(outliers==0));
    stdev = nanstd(data(outliers==0));
else
    mid=mean(data(outliers==0),'omitnan');
    stdev=std(data(outliers==0),'omitnan');
end
data_min = min(data);
data_max = max(data);

[maxFoM, maxFoMat] = max(fom.*(1-outliers));
[maxresult, maxresultat] = max(data.*(1-outliers));



plot_min = data_min;
if handles.trimnegative
    plot_min = mid-handles.varno*stdev;
end

plot_max = data_max;
if handles.trimpositive
    plot_max = mid+handles.varno*stdev;
end

range= plot_max-plot_min;

plot_min = -0.05*range+(plot_min * handles.zoomfrac + data_min * (1-handles.zoomfrac));
plot_max = 0.05*range+(plot_max * handles.zoomfrac + data_max * (1-handles.zoomfrac));


scatter(handles.xydata(handles.xnum,outliers==0),handles.xydata(handles.ynum,outliers==0),18,data(outliers==0),'parent', handles.scatterplot,'filled');
hold(handles.scatterplot, 'on');

if handles.plotoutliers
    scatter(handles.xydata(handles.xnum,outliers==1),handles.xydata(handles.ynum,outliers==1),'b*','parent', handles.scatterplot);
end

if handles.nantoggle
    nanpoint = isnan(data);
    scatter(handles.xydata(handles.xnum,nanpoint),handles.xydata(handles.ynum,nanpoint),'r*','parent', handles.scatterplot);
end

if handles.highlighttopresult
   scatter(handles.xydata(handles.xnum,maxresultat), handles.xydata(handles.ynum,maxresultat),80,'g','filled','parent', handles.scatterplot);
end

if handles.highlighttopFoM
   scatter(handles.xydata(handles.xnum,maxFoMat), handles.xydata(handles.ynum,maxFoMat),80,'r','filled','parent', handles.scatterplot);
end

hold(handles.scatterplot, 'off');



if(plot_min~= plot_max && ~isnan(plot_min) && ~isnan(plot_max))
    set(handles.scatterplot,'clim', [plot_min,plot_max]);
end
set(get(handles.scatterplot,'XLabel'),'string', (handles.xname));
set(get(handles.scatterplot,'YLabel'),'string', (handles.yname));


scatter(handles.xydata(handles.xnum,outliers==0),data(outliers==0), 10, fom(outliers==0),'parent', handles.xplot,'filled');
hold(handles.xplot, 'on');

if handles.plotoutliers
    scatter(handles.xydata(handles.xnum,outliers==1),data(outliers==1),'b*','parent', handles.xplot);
end

if handles.highlighttopresult
   scatter(handles.xydata(handles.xnum,maxresultat), data(maxresultat),80,'g','filled','parent', handles.xplot);
end

if handles.highlighttopFoM
   scatter(handles.xydata(handles.xnum,maxFoMat), data(maxFoMat),80,'r','filled','parent', handles.xplot);
end

hold(handles.xplot, 'off');


if(plot_min~= plot_max && ~isnan(plot_min) && ~isnan(plot_max))
    set(handles.xplot,'ylim', [plot_min,plot_max]);
end
set(get(handles.xplot,'XLabel'),'string', (handles.xname));
set(get(handles.xplot,'YLabel'),'string', (handles.zname));

scatter(handles.xydata(handles.ynum,(outliers==0)),data(outliers==0), 10, fom(outliers==0),'parent', handles.yplot,'filled');

hold(handles.yplot, 'on');
if handles.plotoutliers
    scatter(handles.xydata(handles.ynum,outliers==1),data(outliers==1),'b*','parent', handles.yplot);
end

if handles.highlighttopresult
   scatter(handles.xydata(handles.ynum,maxresultat), data(maxresultat),80,'g','filled','parent', handles.yplot);
end

if handles.highlighttopFoM
   scatter(handles.xydata(handles.ynum,maxFoMat), data(maxFoMat),80,'r','filled','parent', handles.yplot);
end

    hold(handles.yplot, 'off');

if(plot_min~= plot_max && ~isnan(plot_min) && ~isnan(plot_max))
    set(handles.yplot,'ylim', [plot_min,plot_max]);
end
set(get(handles.yplot,'XLabel'),'string', (handles.yname));
set(get(handles.yplot,'YLabel'),'string', (handles.zname));


% --- Executes on button press in datatrimlow.
function datatrimlow_Callback(hObject, eventdata, handles)
% hObject    handle to datatrimlow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of datatrimlow

handles.trimnegative = get(hObject,'Value');
guidata(hObject, handles);
drawgraph(handles)

function stdno_Callback(hObject, eventdata, handles)
% hObject    handle to stdno (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stdno as text
%        str2double(get(hObject,'String')) returns contents of stdno as a double

try
    handles.varno = str2double(get(hObject,'String'));
    if handles.varno<0
        error('Must be positive');
    end
catch
    set(hObject,'String','1');
    handles.varno = 1;
end

guidata(hObject, handles);
drawgraph(handles)



% --- Executes during object creation, after setting all properties.
function stdno_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stdno (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in datatrimhigh.
function datatrimhigh_Callback(hObject, eventdata, handles)
% hObject    handle to datatrimhigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of datatrimhigh

handles.trimpositive = get(hObject,'Value');
guidata(hObject, handles);
drawgraph(handles)


% --- Executes on button press in plotoutliers.
function plotoutliers_Callback(hObject, eventdata, handles)
% hObject    handle to plotoutliers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotoutliers
handles.plotoutliers = get(hObject,'Value');
guidata(hObject, handles);
drawgraph(handles)


% --- Executes on button press in TopFoM.
function TopFoM_Callback(hObject, eventdata, handles)
% hObject    handle to TopFoM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.highlighttopFoM = get(hObject,'Value');
guidata(hObject, handles);
drawgraph(handles)

% --- Executes on button press in TopResult.
function TopResult_Callback(hObject, eventdata, handles)
% hObject    handle to TopResult (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.highlighttopresult = get(hObject,'Value');
guidata(hObject, handles);
drawgraph(handles)


% --- Executes on button press in NaNBox.
function NaNBox_Callback(hObject, eventdata, handles)
% hObject    handle to NaNBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.nantoggle = get(hObject,'Value');
guidata(hObject, handles);
drawgraph(handles)