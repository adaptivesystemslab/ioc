function str = arrayToStr(array, delim)

if nargin == 1
    delim = '-';
end

str = '';

for i = 1:length(array)
    str = [str delim num2str(array(i))];
end
str = str(2:end);
end
