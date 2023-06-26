function saveInput(time,input,temp,SOC,filename)

% template = '../XLSX_INPUTS/Template_input.xlsx';
% fileName = sprintf('../XLSX_INPUTS/%s.xlsx',filename);
template = 'XLSX_INPUTS/Template_input.xlsx';
fileName = sprintf('XLSX_INPUTS/%s.xlsx',filename);

% Prepare new workbook, determine parameter row positions in workbook.
copyfile(template,fileName)

rangeTime  = sprintf('A4:A%d',4+length(input));
rangeInput = sprintf('B4:B%d',4+length(input));
rangeTemp  = sprintf('C4:C%d',4+length(input));
cellSOC    = 'E4';
% writecell(double(time),fileName,'Sheet','Input','Range','A4:A105')

writetable(table(time),fileName, ...
          'Sheet','Input', ...
          'Range',rangeTime, ...
          'AutoFitWidth',false, ...
          'WriteVariableNames',0);

writetable(table(input),fileName, ...
          'Sheet','Input', ...
          'Range',rangeInput, ...
          'AutoFitWidth',false, ...
          'WriteVariableNames',0);

writetable(table(temp),fileName, ...
          'Sheet','Input', ...
          'Range',rangeTemp, ...
          'AutoFitWidth',false, ...
          'WriteVariableNames',0);

writetable(table(SOC),fileName, ...
          'Sheet','Input', ...
          'Range',cellSOC, ...
          'AutoFitWidth',false, ...
          'WriteVariableNames',0);

end