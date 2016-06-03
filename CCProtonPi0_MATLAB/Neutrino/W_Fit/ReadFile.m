function [ x,y ] = ReadFile( file_name )

data = load(file_name);

x = data(:,1);
y = data(:,2);
end

