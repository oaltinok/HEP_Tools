function [ p ] = get_fit_pars( source_file )

[x, y] = get_values(source_file);

p = polyfit(x,y,1);

create_plot(x,y,p);

end

