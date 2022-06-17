 function warn(varargin)
%function warn(varargin)
%
% print concise warning message

[name line] = caller_name(1);
fprintf(['Warn: %s %d: ' varargin{1}], name, line, varargin{2:end})
