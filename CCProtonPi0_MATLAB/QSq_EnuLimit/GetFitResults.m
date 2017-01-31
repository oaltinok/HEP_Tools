function [  ] = GetFitResults(var_name)

fname = strcat(var_name,'.txt');
data = load(fname);
x = data(:,1);
y = data(:,2);

f = fit(x,y,'exp1')
confint(f,0.6827)

end

