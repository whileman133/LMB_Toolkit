function isModel = isCellModel(inModel)
%GETCELLMODELTYPE Determine if the input is a cell model.

throwError = false;
inModelType = getCellModelType(inModel,throwError);
if isempty(inModelType)
    isModel = false;
else
    isModel = true;
end

end