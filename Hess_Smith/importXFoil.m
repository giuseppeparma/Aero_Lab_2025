function importFile=importXFoil(opts,filename)
%
% This function reads XFoil output files and imports the data on Matlab
% This function allows to choose which file to import thanks to the
% switch-case function (according to the output format, different files
% must be read in different ways)
% 
% INPUT
% Name           Type           Size
% opts           char           -
% filename       char           -
%
% OPTIONS
% 'profile'      reads txt file with x,y coordinates describing a profile      
% 'CP'           reads txt file  with computed pressure coeffiicient values 
%                (both viscous and inviscid mode))
% 'CF'           reads txt file with computed friction coefficient values
% 'Coeffs'       reads txt file with a table of commputed aerodynamic
%                coefficients and coordinates of transition points
%
% OUTPUT
% Name           Type           Size
% importFile     struct         1x1
%


switch opts

    case 'profile'
    % Definition of rows to read
    startRow = 2;
    endRow = inf;
    
    % Line format spec
    formatSpec = '%16f%f%[^\n\r]';
    
    % File opening
    fileID = fopen(filename,'r');
    
    % File reading
    dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for block=2:length(startRow)
        frewind(fileID);
        dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
        for col=1:length(dataArray)
            dataArray{col} = [dataArray{col};dataArrayBlock{col}];
        end
    end
    
    % File closing
    fclose(fileID);
    
    importFile = table(dataArray{1:end-1}, 'VariableNames', {'x','y'});



    case 'CP'
        % Definition of rows to read 
    startRow = 2;
    endRow = inf;
    
    % Line format spec
    formatSpec = '%16f%f%[^\n\r]';
    
    % File opening
    fileID = fopen(filename,'r');
    
    % File reading
    dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for block=2:length(startRow)
        frewind(fileID);
        dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
        for col=1:length(dataArray)
            dataArray{col} = [dataArray{col};dataArrayBlock{col}];
        end
    end
    
    % File closing
    fclose(fileID);

    importFile = table(dataArray{1:end-1}, 'VariableNames', {'x','Cp'});



    case 'CF'
    % Definition of rows to read
    startRow = 8;
    endRow = inf;
    
    
    % Line format spec
    formatSpec = '%16f%f%[^\n\r]';
    
    % File opening
    fileID = fopen(filename,'r');
    
    % File reading
    dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for block=2:length(startRow)
        frewind(fileID);
        dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
        for col=1:length(dataArray)
            dataArray{col} = [dataArray{col};dataArrayBlock{col}];
        end
    end
    
    % File closing
    fclose(fileID);

    importFile = table(dataArray{1:end-1}, 'VariableNames', {'x','Cf'});
    
    
    
    case 'Coeffs'
    % Definition of rows to read
    startRow = 13;
    endRow=Inf;
    
    % LLLine format spec
    formatSpec = '%16f%f%[^\n\r]';

    % Open the text file.
    fileID = fopen(filename,'r');

    % Skip lines up to startRow-1: this process perrforms the fgetl the number
    % of times required to reach the actual startRow

    i=1;
    for k = 1:startRow-1
        if feof(fileID), break; end %exits if the end of the file is reached (i.e. if the file is shorter than the prescribed startRow number)
        fgetl(fileID); % moves to the next line
    end

    % Now read one row at a time and process
    while ~feof(fileID)
        % Scan one row at a a time
        C = textscan(fileID, formatSpec, 1, ...
                 'Delimiter','', 'WhiteSpace',' \t', ...
                 'EndOfLine','\n', 'ReturnOnError', false, 'TextType','string');
        % If textscan returned empty (EOF or bad parse), break
        if all(cellfun(@(x) isempty(x), C)), break; end
        if i>endRow, break; end
        % Extract and process row values (the third column mmust be processed to have each of its values placed in a separate cell)
        a = C{1}; b = C{2}; 
        s = C{3};
        
        % Divide third cell to have an arrray with one numeric value in each cell
        nums = str2num(s);           
        out = num2cell(nums);        
        
        % Build the output in the right form
        data(i,:)=[a, b, out];
        i = i + 1; % Increment the row index for the next iteration
    
    end

    % Close the text file.
    fclose(fileID);

    % Convert cell array into mattrix
    data=cell2mat(data);

    % Create output variable
    importFile = table(data(:,1),data(:,2),data(:,3),data(:,4),data(:,5),data(:,6),data(:,7),data(:,8),data(:,9), 'VariableNames', {'alpha','cl','cd','cdp','cm','top_xtr','bot_xtr','top_itr','bot_itr'});

end