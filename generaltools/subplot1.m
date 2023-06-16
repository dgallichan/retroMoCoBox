function hAxes = subplot1(M,N,varargin);
%-------------------------------------------------------------------------
% subplot1 function         An mproved subplot function
% Input  : - If more than one input argumenst are given,
%            then the first parameter is the number of rows.
%            If single input argument is given, then this is the
%            subplot-number for which to set focus.
%            This could a scalar or two element vector (I,J).
%          - Number of columns.
%          * variable number of parameters
%            (in pairs: ...,Keywoard, Value,...)
%           - 'Min'    : X, Y lower position of lowest subplot,
%                        default is [0.10 0.10].
%           - 'Max'    : X, Y largest position of highest subplot,
%                        default is [0.95 0.95].
%           - 'Gap'    : X,Y gaps between subplots,
%                        default is [0.01 0.01].
%           - 'XTickL' : x ticks labels option,
%                        'Margin' : plot only XTickLabels in the
%                                   subplot of the lowest  row (default).
%                        'All'    : plot XTickLabels in all subplots.
%                        'None'   : don't plot XTickLabels in subplots.
%           - 'YTickL' : y ticks labels option,
%                        'Margin' : plot only YTickLabels in the
%                                   subplot of the lowest  row (defailt).
%                        'All'    : plot YTickLabels in all subplots.
%                        'None'   : don't plot YTickLabels in subplots.
%           -  'FontS'  : axis font size, default is 10.
%             'XScale' : scale of x axis:
%                        'linear', default.
%                        'log'
%           -  'YScale' : scale of y axis:
%                        'linear', default.
%                        'log'
% Example: subplot1(2,2,'Gap',[0.02 0.02]);
%          subplot1(2,3,'Gap',[0.02 0.02],'XTickL','None','YTickL','All','FontS',16);
% See also : subplot1c.m
% Tested : Matlab 5.3
%     By : Eran O. Ofek           June 2002
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%-------------------------------------------------------------------------

MinDef      = [0.10 0.10];
MaxDef      = [0.95 0.95];
GapDef      = [0.02 0.01];
XTickLDef   = 'Margin';  
YTickLDef   = 'Margin';  
FontSDef    = 10;
XScaleDef   = 'linear';
YScaleDef   = 'linear';

figHandleDef = [];

% set default parameters
Min    = MinDef;
Max    = MaxDef;
Gap    = GapDef;
XTickL = XTickLDef;
YTickL = YTickLDef;
FontS  = FontSDef;
XScale = XScaleDef;
YScale = YScaleDef;

figHandle = figHandleDef;

MoveFoc = 0;
if (nargin==1),
   %--- move focus to subplot # ---
   MoveFoc = 1;
elseif (nargin==2),
   % do nothing
elseif (nargin>2),
   Narg = length(varargin);
   if (0.5.*Narg==floor(0.5.*Narg)),

      for I=1:2:Narg-1,
         switch varargin{I},
          case 'Min'
 	     Min = varargin{I+1};
          case 'Max'
 	     Max = varargin{I+1};
          case 'Gap'
 	     Gap = varargin{I+1};
          case 'XTickL'
 	     XTickL = varargin{I+1};
          case 'YTickL'
 	     YTickL = varargin{I+1};
          case 'FontS'
 	     FontS = varargin{I+1};
          case 'XScale'
 	     XScale = varargin{I+1};
          case 'YScale'
 	     YScale = varargin{I+1};
          case 'figHandle' % danielg
         figHandle = varargin{I+1};                      
          otherwise
	     error('Unknown keyword');
         end
      end
   else
      error('Optional arguments should given as keyword, value');
   end
else
   error('Illegal number of input arguments');
end




switch MoveFoc
    case 1
        % danielg:
        if isempty(figHandle)
            figHandle = gcf;
        end
        
        %--- move focus to subplot # ---
        H    = get(figHandle,'Children');
        %     Ptot = length(H);
        iReject = [];
        for iP = 1:length(H) %% danielg: try to fix for HG2
            if strcmp(get(H(iP),'tag'),'Colorbar')
                iReject = [iReject iP];
            end
        end
        H(iReject) = [];
        
        Ptot = sum(strcmp(get(H,'type'),'axes'));
        if (length(M)==1),
            M    = Ptot - M + 1;
        elseif (length(M)==2),
            %--- check for subplot size ---
            Pos1  = get(H(1),'Position');
            Pos1x = round(Pos1(1)*100)/100;
            for Icheck=2:1:Ptot,
                PosN  = get(H(Icheck),'Position');
                PosNx = round(PosN(1)*100)/100;
                if (PosNx==Pos1x),
                    NumberOfCol = Icheck - 1;
                    break;
                end
            end
            NumberOfRow = Ptot./NumberOfCol;
            
            Row = M(1);
            Col = M(2);
            
            M   = (Row-1).*NumberOfCol + Col;
            M    = Ptot - M + 1;
        else
            error('Unknown option, undefined subplot index');
        end

        if nargout > 0  hAxes = H(M);  else set(figHandle,'CurrentAxes',H(M)); end
 

    case 0
        %--- open subplots ---
        
        Xmin   = Min(1);
        Ymin   = Min(2);
        Xmax   = Max(1);
        Ymax   = Max(2);
        Xgap   = Gap(1);
        Ygap   = Gap(2);
        
        
        Xsize  = (Xmax - Xmin)./N;
        Ysize  = (Ymax - Ymin)./M;
        
        Xbox   = Xsize - Xgap;
        Ybox   = Ysize - Ygap;
        
        
        Ptot = M.*N;
        
        if isempty(figHandle)
            Hgcf = gcf;
        else
            Hgcf = figHandle;
        end
        clf;
        fig(Hgcf);
        for Pi=1:1:Ptot,
            Row = ceil(Pi./N);
            Col = Pi - (Row - 1)*N;
            
            Xstart = Xmin + Xsize.*(Col - 1);
            Ystart = Ymax - Ysize.*Row;
            
            if nargout > 0
                hAxes(Row,Col) = axes('position',[Xstart,Ystart,Xbox,Ybox]);
            else
                axes('position',[Xstart,Ystart,Xbox,Ybox]);
            end
            
            %       subplot(M,N,Pi);
            %       hold on;
            %set(gca,'position',[Xstart,Ystart,Xbox,Ybox]);
            set(gca,'FontSize',FontS);
            
            set(gca,'Tag',num2str(Pi));
            
            box on;
            hold on;
            
            switch XTickL
                case 'Margin'
                    if (Row~=M),
                        %--- erase XTickLabel ---
                        set(gca,'XTickLabel',[]);
                    end
                case 'All'
                    % do nothing
                case 'None'
                    set(gca,'XTickLabel',[]);
                otherwise
                    error('Unknown XTickL option');
            end
            
            switch YTickL
                case 'Margin'
                    if (Col~=1),
                        %--- erase YTickLabel ---
                        set(gca,'YTickLabel',[]);
                    end
                case 'All'
                    % do nothing
                case 'None'
                    set(gca,'YTickLabel',[]);
                otherwise
                    error('Unknown XTickL option');
            end
            
            switch XScale
                case 'linear'
                    set(gca,'XScale','linear');
                case 'log'
                    set(gca,'XScale','log');
                otherwise
                    error('Unknown XScale option');
            end
            
            switch YScale
                case 'linear'
                    set(gca,'YScale','linear');
                case 'log'
                    set(gca,'YScale','log');
                otherwise
                    error('Unknown YScale option');
            end
            
        end
        
    otherwise
        error('Unknown MoveFoc option');
end
