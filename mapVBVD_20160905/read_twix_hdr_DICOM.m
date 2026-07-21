function prot = read_twix_hdr_DICOM(DICOMfile)
% modified version of Philip Ehses' read_twix_hdr.m to run on DICOM files
% and extract the Siemens parameters from the Private_1029_1020 tag
%
% It is not clear if the XML part of this header is actually useful - but
% the ASCCONV section clearly contains all the 'regular' Siemens parameters
% used to define the sequence. Extracting the XML doesn't work easily with
% the sam fucntionality as read_twix_hdr.m, so for now this just gets the
% 'ASCCONV' part.
%
     

   dcminfo = dicominfo(DICOMfile);
   str=char(dcminfo.Private_0029_1020(:))';
   
   [ascconv,xprot_noUse] = regexp(str,'### ASCCONV BEGIN[^\n]*\n(.*)\s### ASCCONV END ###','tokens','split');
  
   if ~isempty(ascconv)
       ascconv = ascconv{:}{:};
       prot = parse_ascconv(ascconv);
   else
       prot = [];
   end
   
   %%% None of the experiments below seemed to be able to extract anything
   %%% useful...
%    iStartX = strfind(str,'<XProtocol>');
%    iStartAsc = strfind(str,'### ASCCONV BEGIN');
%    
%    prot.xprot1 = parse_xprot(str(iStartX(1):iStartX(2)-1));
%    prot.xprot2 = parse_xprot(str(iStartX(2):iStartX(3)-1));
%    prot.xprot3 = parse_xprot(str(iStartX(3):iStartAsc-1));

end
   
   
       
    


function xprot = parse_xprot(buffer)
    xprot = [];
    tokens = regexp(buffer, '<Param(?:Bool|Long|String)\."(\w+)">\s*{([^}]*)','tokens');
    tokens = [tokens, regexp(buffer, '<ParamDouble\."(\w+)">\s*{\s*(<Precision>\s*[0-9]*)?\s*([^}]*)','tokens')];
    for m=1:numel(tokens)
        name         = char(tokens{m}(1));
        % field name has to start with letter
        if (~isletter(name(1)))
            name = strcat('x', name);
        end

        value = char(strtrim(regexprep(tokens{m}(end), '("*)|( *<\w*> *[^\n]*)', '')));
        value = regexprep(value, '\s*', ' ');

        try %#ok<TRYNC>
            value = eval(['[' value ']']);  % inlined str2num()
        end

        xprot.(name) = value;
    end
end

function mrprot = parse_ascconv(buffer)  
    mrprot = [];    
    % [mv] was: vararray = regexp(buffer,'(?<name>\S*)\s*=\s(?<value>\S*)','names');
    vararray = regexp(buffer,'(?<name>\S*)\s*=\s*(?<value>\S*)','names');
    
    for var=vararray

        try
            value = eval(['[' var.value ']']);  % inlined str2num()
        catch
            value = var.value;
        end
        
        % now split array name and index (if present)
        v = regexp(var.name,'(?<name>\w*)\[(?<ix>[0-9]*)\]|(?<name>\w*)','names');

        cnt = 0;
        tmp = cell(2,numel(v));

        breaked = false;
        for k=1:numel(v)
            if ~isletter(v(k).name(1))
                breaked = true;
                break;
            end
            cnt = cnt+1;
            tmp{1,cnt} = '.';
            tmp{2,cnt} = v(k).name;

            if ~isempty(v(k).ix)
                cnt = cnt+1;
                tmp{1,cnt} = '{}';
                tmp{2,cnt}{1} = 1 + str2double(v(k).ix);
            end
        end
        if ~breaked && ~isempty(tmp)
            S = substruct(tmp{:});
            mrprot = subsasgn(mrprot,S,value);
        end
    end 
end

