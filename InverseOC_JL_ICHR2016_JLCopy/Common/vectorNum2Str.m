function str = vectorNum2Str(vector, lengthVecString)
  str = ['['];
  
  for i = 1:length(vector)
      str = [str num2str(vector(i), lengthVecString) ','];
  end
  
  str = [str(1:end-1) ']'];
end