function [x,y]=createXfoil(opts,naca,NPanels,Chord,alpha,Re,iter)

% createXfoil
% This function creates the input file for XFoil and executes it 
% This function allows to choose what to do on Xfoil thanks to the
% switch-case function
%
% INPUT
% Name           Type           Size
% opts           char           -
% naca           char           -
% NPanels        scalar         1x1
% Chord          scalar         1x1
% alpha          scalar         1x1
%                vector         3x1
%
% if alpha is a vector, it represents a sequence of values and implies an analysis for each of them;
% syntax: alpha=(firstValue,lastValue,increment)
% firstValue    first values of alpha in the sequence
% lastValue     last value of alpha  in the sequence
% increment     increment to be used to go from one value to the next in the sequence
%
% Re             scalar         1x1
% iter           scalar         1x1
%
% OPTIONS
% 'profile'      creates the profile from XFoil
% syntax         [x,y]=createXfoil('profile',naca,NPanels,Chord)
% 
% 'CP'           creates the viscous pressure coefficient from XFoil
% syntax         [x,y]=createXfoil('CP',naca,NPanels,Chord,alpha)
% 
% 'CPvisc'       creates the pressure coefficient from XFoil
% syntax         [x,y]=createXfoil('CPvisc',naca,NPanels,Chord,alpha,Re,iter)
% 
% 'CF'           creates the friction coefficient from XFoil
% synatx         [x,y]=createXfoil('CF',naca,NPanels,Chord,alpha,Re,iter)
% 
% 'Coeffs'       creates the coeffficients related to prescibed alpha values 
% syntax         [x]=createXfoil('coeffs',naca,NPanels,Chord,alpha)
% 
% 'Coeffsvisc'   creates the coeffficients in viscous mode related to prescibed alpha values 
% syntax         [x]=createXfoil('coeffs',naca,NPanels,Chord,alpha,Re)
% 
% 
% OUTPUT
% Name           Type           Size
% x              vector         nRowsxnCols
%                struct         nRowsx9
% y              vector         nRowsxnCols
%
% the output dimensions vary according to the option requested when calling
% the function

switch opts

    case 'profile'
        % Creation of Xfoil  input file
        fileID = fopen('XFoilInput.txt','w');
        fprintf(fileID, ['naca' ' '  naca, '\n\n']); 
        fprintf(fileID,'pane\n\n'); 

        fprintf(fileID,'gdes\n');
        fprintf(fileID,'tgap 0 0 \n');
    
        fprintf(fileID,'exec \n\n\n');

        fprintf(fileID,'ppar\n');
        fprintf(fileID, ['n ' ' ' num2str(NPanels+1)  '\n\n\n']);

        filename = strcat('NACA_', naca, '.txt');

        fprintf(fileID, ['save ' ' ' filename '\n\n']);
        fprintf(fileID,'y\n\n');
        fprintf(fileID,'quit \n\n');
        fclose(fileID);

        Str2Exec = strcat("xfoil < XFoilInput.txt > /dev/null 2>&1");
        system(Str2Exec);
        
        % Import data from Xfoill ooutput file
        profile = importXFoil('profile',filename);
    
        % Vectors inversion
        % output dimensions: nRows=NPanels, nCols=1;
        x = flipud(profile.x);
        y = flipud(profile.y);
      
        x = x.*Chord;
        y = y.*Chord;


    case 'CP'
        % Creation of Xfoil input file
        fileID = fopen('XFoilCPInput.txt','w');
        fprintf(fileID, ['naca' ' '  naca, '\n\n']); 
        fprintf(fileID,'pane\n\n'); 

        fprintf(fileID,'gdes\n');
        fprintf(fileID,'tgap 0 0 \n');
        
        fprintf(fileID,'exec \n\n\n');

        fprintf(fileID,'ppar\n');
        fprintf(fileID, ['n ' ' ' num2str(NPanels+1)  '\n\n\n']);

        fprintf(fileID, 'oper\n');
        fprintf(fileID, 'alfa\n');
        fprintf(fileID, [' ' num2str(alpha), '\n']);
        
        filename = strcat('NACA_CP_', naca, '.txt');
        
        fprintf(fileID, ['cpwr ' ' ' filename '\n\n']);
        fprintf(fileID,'y\n\n');
        fprintf(fileID,'quit \n\n');
        fclose(fileID);

        Str2Exec = strcat("xfoil < XFoilCPInput.txt > /dev/null 2>&1");
        system(Str2Exec);
        
        % Import data from Xffoil output file
        profile = importXFoil('CP',filename);
       
        % Vectors inversion
        % output dimensions: nRows=NPanels, nCols=1;
        x_Cp = flipud(profile.x);
        Cp = flipud(profile.Cp);

        x = x_Cp.*Chord;
        y = Cp.*Chord;

    case 'CPvisc'
        % Creation of Xfoil input file
        fileID = fopen('XFoilCPviscInput.txt','w');
        fprintf(fileID, ['naca' ' '  naca, '\n\n']); 
        fprintf(fileID,'pane\n\n'); 

        fprintf(fileID,'gdes\n');
        fprintf(fileID,'tgap 0 0 \n');
        
        fprintf(fileID,'exec \n\n\n');

        fprintf(fileID,'ppar\n');
        fprintf(fileID, ['n ' ' ' num2str(NPanels+1)  '\n\n\n']);

        fprintf(fileID, 'oper\n');
        fprintf(fileID, 'visc\n');
        fprintf(fileID, [' ' num2str(Re) '\n']);
        fprintf(fileID, ['iter ' num2str(iter) '\n']);
        fprintf(fileID, 'alfa\n');
        fprintf(fileID, [' ' num2str(alpha), '\n']);
        
        filename = strcat('NACA_CPvisc_', naca, '.txt');
        
        fprintf(fileID, ['cpwr ' ' ' filename '\n\n']);
        fprintf(fileID,'y\n\n');
        fprintf(fileID,'quit \n\n');
        fclose(fileID);

        Str2Exec = strcat("xfoil < XFoilCPviscInput.txt > /dev/null 2>&1");
        system(Str2Exec);

        % Import data from Xfoil output file
        profile = importXFoil('CP',filename);
       
        % Vectors inversion
        % output dimensions: nRows=NPanels, nCols=1;
        x_Cp = flipud(profile.x);
        Cp = flipud(profile.Cp);

        x = x_Cp.*Chord;
        y = Cp.*Chord;

    case 'CF'
            
        filename = strcat('NACAcf_', naca, '.txt');
        % Delete file if already existing to avoid overwriting errors
        if isfile(filename)
            delete(filename);
        end
        
        % Creation of Xfoil input file    
        fileID = fopen('XFoilCfInput.txt','w');
        fprintf(fileID, ['naca' ' '  naca, '\n\n']); 
        fprintf(fileID,'pane\n\n'); 

        fprintf(fileID,'gdes\n');
        fprintf(fileID,'tgap 0 0 \n');
        
        fprintf(fileID,'exec \n\n\n');

        fprintf(fileID,'ppar\n');
        fprintf(fileID, ['n ' ' ' num2str(NPanels+1)  '\n\n\n']);

        fprintf(fileID, 'oper\n');
        fprintf(fileID, 'visc\n');
        fprintf(fileID, [' ' num2str(Re) '\n']);
        fprintf(fileID, ['iter ' ' ' num2str(iter) '\n']);

        fprintf(fileID, 'alfa\n');
        fprintf(fileID, [' ' num2str(alpha), '\n']);
        
        fprintf(fileID, 'vplo\n');
        fprintf(fileID, 'cf\n');

        fprintf(fileID, ['dump ' ' ' filename '\n\n']);
        
        
        fprintf(fileID,'y\n\n');
        fprintf(fileID,'quit \n\n');
        fclose(fileID);

        Str2Exec = strcat("xfoil < XFoilCfInput.txt > /dev/null 2>&1");
        system(Str2Exec);

        %  Import data from Xfoil output file
        File = importXFoil('CF',filename);
       
        % Process imported data to obtain output matrices
        % The structure of the Xfoil output file requires post-processing 
        % to divide the coefficient values on the suction side from the 
        % ones on the pressure side
  
        tic=[];
        for i=1:207
            test=File.x(i);
            if test<=1 
                importFile(i,1)=File.x(i);
                importFile(i,2)=File.Cf(i);
            else 
                tic=[tic,i];
            end
        end

        for j=1:tic(1)-1
            importTop(j,1)=importFile(j,1);
            importTop(j,2)=importFile(j,2);
        end

        l=1;
        toc=length(tic);
        k1=tic(toc/2)+1;
        for k=k1:length(importFile(:,1))
            importBot(l,1)=importFile(k,1);
            importBot(l,2)=importFile(k,2);
            l=l+1;
        end        
        
        % output dimensions: nRowsx2 matrices
        % nRows are related to the number of values imported from Xfoil
        x = importTop;
        y = importBot;

    case 'Coeffs'
        filename = strcat('NACAcoefficients_', naca, '.txt');
        % Delete file if already existing to  avoid overwriting errors
        if isfile(filename)
            delete(filename);
        end
        
       
        fileID = fopen('XFoilCoefficientsInput.txt','w');
        fprintf(fileID, ['naca' ' '  naca, '\n\n']); 
        fprintf(fileID,'pane\n\n'); 

        fprintf(fileID,'gdes\n');
        fprintf(fileID,'tgap 0 0 \n');
        
        fprintf(fileID,'exec \n\n\n');

        fprintf(fileID,'ppar\n');
        fprintf(fileID, ['n ' ' ' num2str(NPanels+1)  '\n\n\n']);
        
       
        fprintf(fileID, 'oper\n');
        fprintf(fileID, 'seqp\n');
        fprintf(fileID, 'pacc\n');
        fprintf(fileID, [' ' filename '\n\n']);

        % Distinguish between case of single alpha value and alpha sequence
        if isscalar(alpha)
            fprintf(fileID, ['alfa ' ' ' num2str(alpha) '\n\n\n']);
        
        else
            fprintf(fileID, 'aseq \n');
            fprintf(fileID, [' ' num2str(alpha(1)) '\n']);
            fprintf(fileID, [' ' num2str(alpha(2)) '\n']);
            fprintf(fileID, [' ' num2str(alpha(3)) '\n']);
        end
        

        fprintf(fileID,'y\n\n');
        fprintf(fileID,'quit \n\n');
        fclose(fileID);

        Str2Exec = strcat("xfoil < XFoilCoefficientsInput.txt > /dev/null 2>&1");

        system(Str2Exec);

        % Import data from Xfoil output file 
        Corpo = importXFoil('Coeffs',filename);
       
        % output dimensions: struct of 9 fields
        % nRows=number of alpha values in alpha sequence and nCols=9
        x = Corpo;
        y = 0;


    case 'Coeffsvisc'

        filename = strcat('NACAviscCoefficients_', naca, '.txt');
        % Delete file if already existing to  avoid overwriting errors
        if isfile(filename)
            delete(filename);
        end

        

        fileID = fopen('XFoilViscCoefficientsInput.txt','w');
        fprintf(fileID, ['naca' ' '  naca, '\n\n']); 
        fprintf(fileID,'pane\n\n'); 

        fprintf(fileID,'gdes\n');
        fprintf(fileID,'tgap 0 0 \n');

        fprintf(fileID,'exec \n\n\n');

        fprintf(fileID,'ppar\n');
        fprintf(fileID, ['n ' ' ' num2str(NPanels+1)  '\n\n\n']);


        fprintf(fileID, 'oper\n');
        fprintf(fileID, 'visc\n');
        fprintf(fileID, [' ' num2str(Re) '\n']);
        fprintf(fileID, 'iter\n');
        fprintf(fileID, [' ' num2str(iter) '\n']);
        fprintf(fileID, 'seqp\n');
        fprintf(fileID, 'pacc\n');
        fprintf(fileID, [' ' filename '\n\n']);

        % Distinguish between case of single alpha value and alpha sequence
        if isscalar(alpha)
            fprintf(fileID, ['alfa ' ' ' num2str(alpha) '\n\n\n']);

        else
            fprintf(fileID, 'aseq \n');
            fprintf(fileID, [' ' num2str(alpha(1)) '\n']);
            fprintf(fileID, [' ' num2str(alpha(2)) '\n']);
            fprintf(fileID, [' ' num2str(alpha(3)) '\n']);
        end


        fprintf(fileID,'y\n\n');
        fprintf(fileID,'quit \n\n');
        fclose(fileID);

        Str2Exec = strcat("xfoil < XFoilViscCoefficientsInput.txt > /dev/null 2>&1");

        system(Str2Exec);

        % Import data from Xfoil output file 
        Corpo = importXFoil('Coeffs',filename);

        % output dimensions: struct of 9 fields
        % nRows=number of alpha values in alpha sequence and nCols=9
        x=Corpo;

        y = 0;

end

end
