function [yout,x] = ow_imhist(varargin)
%IMHIST Display histogram of image data.
%   IMHIST(I) displays a histogram for the intensity image I whose number of
%   bins are specified by the image type.  If I is a grayscale image, IMHIST
%   uses 256 bins as a default value. If I is a binary image, IMHIST uses
%   only 2 bins.
%
%   IMHIST(I,N) displays a histogram with N bins for the intensity image I
%   above a grayscale colorbar of length N.  If I is a binary image then N
%   can only be 2.
%
%   IMHIST(X,MAP) displays a histogram for the indexed image X. This
%   histogram shows the distribution of pixel values above a colorbar of the
%   colormap MAP. The colormap must be at least as long as the largest index
%   in X. The histogram has one bin for each entry in the colormap.
%
%   [COUNTS,X] = imhist(...) returns the histogram counts in COUNTS and the
%   bin locations in X so that stem(X,COUNTS) shows the histogram. For
%   indexed images, it returns the histogram counts for each colormap entry;
%   the length of COUNTS is the same as the length of the colormap.
%
%   Class Support
%   -------------  
%   An input intensity image can be uint8, uint16, int16, single, double, or
%   logical. An input indexed image can be uint8, uint16, single, double, or
%   logical.  
%
%   Note
%   ----
%   For intensity images, the N bins of the histogram are each half-open
%   intervals of width A/(N-1).
%  
%   In particular, for intensity images that are not int16, the p-th bin is the
%   half-open interval:
%
%        A*(p-1.5)/(N-1)  <= x  <  A*(p-0.5)/(N-1)
%
%   For int16 intensity images, the p-th bin is the half-open interval:
%  
%        A*(p-1.5)/(N-1) - 32768  <= x  <  A*(p-0.5)/(N-1) - 32768  
%
%   The intensity value is represented by "x". The scale factor A depends on the
%   image class.  A is 1 if the intensity image is double or single; A is 255 if
%   the intensity image is uint8; A is 65535 if the intensity image is uint16 or
%   int16.
%  
%   Example
%   -------
%        I = imread('pout.tif');
%        imhist(I)
%
%   See also HISTEQ, HIST.

%   Copyright 1992-2007 The MathWorks, Inc.
%   $Revision: 5.24.4.8 $ $Date: 2007/07/31 20:02:34 $

[a, n, isScaled, top, map] = parse_inputs(varargin{:});

if islogical(a)
    if (n ~= 2)
        messageId = 'Images:imhist:invalidParameterForLogical';
        message1 = 'N must be set to two for a logical image.'; 
        error(messageId, '%s', message1);
    end
    y(2) = sum(a(:));
    y(1) = numel(a) - y(2);
    y = y';
elseif isa(a,'int16')
    y = imhistc(int16touint16(a), n, isScaled, top); % Call MEX file to do work.
else
    y = imhistc(a, n, isScaled, top); % Call MEX file to do work.
end

range = getrangefromclass(a);
if islogical(a)
    x = range';
else
    % integer or float
    x = linspace(range(1), range(2), n)';
end

if (nargout == 0)
    plot_result(x, y, map, isScaled, class(a), range);
else
    yout = y;
end


%%%
%%% Function plot_result
%%%
function plot_result(x, y, cm, isScaled, classin, range)

n = length(x);
stem(x,y, 'Marker', 'none')
hist_axes = gca;

h_fig = ancestor(hist_axes,'figure');

% Get x/y limits of axes using axis
limits = axis(hist_axes);
if n ~= 1
  limits(1) = min(x);
else
  limits(1) = 0;
end
limits(2) = max(x);
var = sqrt(y'*y/length(y));
limits(4) = 2.5*var;
axis(hist_axes,limits);


% Cache the original axes position so that axes can be repositioned to
% occupy the space used by the colorstripe if nextplot clears the histogram
% axes.
original_axes_pos = get(hist_axes,'Position');

% In GUIDE, default axes units are characters. In order for axes repositiong
% to behave properly, units need to be normalized.
hist_axes_units_old = get(hist_axes,'units');
set(hist_axes,'Units','Normalized');
% Get axis position and make room for color stripe.
pos = get(hist_axes,'pos');
stripe = 0.075;
set(hist_axes,'pos',[pos(1) pos(2)+stripe*pos(4) pos(3) (1-stripe)*pos(4)])
set(hist_axes,'Units',hist_axes_units_old);

set(hist_axes,'xticklabel','')

% Create axis for stripe
stripe_axes = axes('Parent',get(hist_axes,'Parent'),...
                'Position', [pos(1) pos(2) pos(3) stripe*pos(4)]);
				 				 
limits = axis(stripe_axes);

% Create color stripe
if isScaled,
    binInterval = 1/n;
    xdata = [binInterval/2 1-(binInterval/2)];
    limits(1:2) = range;
    switch classin
     case {'uint8','uint16'}
        xdata = range(2)*xdata;
        C = (1:n)/n;
     case 'int16'
        xdata = 65535*xdata - 32768;
        C = (1:n)/n;
     case {'double','single'}
        C = (1:n)/n;
     case 'logical'
        C = [0 1];
     otherwise
        messageId = sprintf('Images:%s:internalError', mfilename);
        error(messageId, 'The input image must be uint8, uint16, %s', ...
            'double, or logical.');    
    end
    
    % image(X,Y,C) where C is the RGB color you specify. 
    image(xdata,[0 1],repmat(C, [1 1 3]),'Parent',stripe_axes);
else
    if length(cm)<=256
        image([1 n],[0 1],1:n,'Parent',stripe_axes); 
        set(h_fig,'Colormap',cm);
        limits(1) = 0.5;
        limits(2) = n + 0.5;
    else
        image([1 n],[0 1],permute(cm, [3 1 2]),'Parent',stripe_axes);
        limits(1) = 0.5;
        limits(2) = n + 0.5;
    end
end

set(stripe_axes,'yticklabel','')
axis(stripe_axes,limits);

% Put a border around the stripe.
line(limits([1 2 2 1 1]),limits([3 3 4 4 3]),...
       'LineStyle','-',...
       'Parent',stripe_axes,...
       'Color',get(stripe_axes,'XColor'));

% Special code for a binary image
if strcmp(classin,'logical')
    % make sure that the stripe's X axis has 0 and 1 as tick marks.
    set(stripe_axes,'XTick',[0 1]);

    % remove unnecessary tick marks from axis showing the histogram
    set(hist_axes,'XTick',0);
    
    % make the histogram lines thicker
    h = get(hist_axes,'children');
    obj = findobj(h,'flat','Color','b');
    lineWidth = 10;
    set(obj,'LineWidth',lineWidth);
end

set(h_fig,'CurrentAxes',hist_axes);

% Tag for testing. 
set(stripe_axes,'tag','colorstripe');

wireHistogramAxesListeners(hist_axes,stripe_axes,original_axes_pos);

% Link the XLim of histogram and color stripe axes together.
% In calls to imhist in a tight loop, the histogram and colorstripe axes
% are destroyed and recreated repetitively. Use linkprop rather than
% linkaxes to link xlimits together to solve deletion timing problems.
h_link = linkprop([hist_axes,stripe_axes],'XLim');
setappdata(stripe_axes,'linkColorStripe',h_link);

%%%
%%% Function wireHistogramAxesListeners
%%%
function wireHistogramAxesListeners(hist_axes,stripe_axes,original_axes_pos)

hist_axes_handle = handle(hist_axes);

% If the histogram axes is deleted, delete the color stripe associated with
% the histogram axes.
lis.histogramAxesDeletedListener = handle.listener(hist_axes_handle,...
                                                'ObjectBeingDestroyed',...
                                                {@removeColorStripeAxes,stripe_axes});

% This is a dummy hg object used to listen for when the histogram axes is cleared.
deleteProxy = handle(text('Parent',hist_axes,...
                          'Visible','Off', ...
                          'Tag','axes cleared proxy',...
                          'HandleVisibility','off'));

% deleteProxy is an invisible text object that is parented to the histogram
% axes.  If the ObjectBeingDestroyed listener fires, the histogram axes has
% been cleared. This listener is triggered by newplot when newplot clears
% the current axes to make way for new hg objects being drawn. This
% listener does NOT fire as a result of the parent axes being deleted.
lis.proxydeleted = handle.listener(deleteProxy,...
                                   'ObjectBeingDestroyed',...
                                   {@histogramAxesCleared,stripe_axes,original_axes_pos});

setappdata(stripe_axes,'ColorStripeListeners',lis);


%%%
%%% Function removeColorStripeAxes
%%%
function removeColorStripeAxes(hDeleteProxy,ed,stripe_axes) %#ok		

if ishandle(stripe_axes)
    delete(stripe_axes);
end
        

%%%
%%% Function histogramAxesCleared
%%%
function histogramAxesCleared(hDeleteProxy,ed,stripe_axes,original_axes_pos)

removeColorStripeAxes(hDeleteProxy,ed,stripe_axes);

h_hist_ax = get(hDeleteProxy,'parent');
set(h_hist_ax,'Position',original_axes_pos);


%%%
%%% Function parse_inputs
%%%
function [a, n, isScaled, top, map] = parse_inputs(varargin)

iptchecknargin(1,2,nargin,mfilename);
a = varargin{1};
iptcheckinput(a, {'double','uint8','logical','uint16','int16','single'}, ...
              {'2d','nonsparse'}, mfilename, ['I or ' 'X'], 1);
n = 256;

if isa(a,'double') || isa(a,'single')
    isScaled = 1;
    top = 1;
    map = []; 
    
elseif isa(a,'uint8')
    isScaled = 1; 
    top = 255;
    map = [];
    
elseif islogical(a)
    n = 2;
    isScaled = 1;
    top = 1;
    map = [];
    
else % int16 or uint16
    isScaled = 1; 
    top = 65535;
    map = [];
end
    
if (nargin ==2)
    if (numel(varargin{2}) == 1)
        % IMHIST(I, N)
        n = varargin{2};
        iptcheckinput(n, {'numeric'}, {'real','positive','integer'}, mfilename, ...
                      'N', 2);
        
    elseif (size(varargin{2},2) == 3)
      if isa(a,'int16')
        eid = sprintf('Images:%s:invalidIndexedImage',mfilename);
        msg1 = 'An indexed image can be uint8, uint16, double, single, or ';
        msg2 = 'logical.';
        error(eid,'%s %s',msg1, msg2);
      end

      % IMHIST(X,MAP) or invalid second argument
      n = size(varargin{2},1);
      isScaled = 0;
      top = n;
      map = varargin{2};
      
    else
        messageId = sprintf('Images:%s:invalidSecondArgument', mfilename);
        message4 = 'Second argument must be a colormap or a positive integer.';
        error(messageId, '%s', message4); 
    end
end
