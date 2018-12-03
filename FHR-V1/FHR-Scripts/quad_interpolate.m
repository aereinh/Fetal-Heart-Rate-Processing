%Given a set of coordinates (x,y) and some input value on x, 
%performs parabolic interpolation of three points surrounding the input
%Returns parabola maximum, x_locs, and coefficients

function [interp_max,x_loc,coeffs] = quad_interpolate(x,y,input)

%Convert complex signals to real ones
if ~isreal(y)
    y = abs(y);
end

%Make sure vectors are column vectors
if isrow(x)
    x = x';
end
if isrow(y)
    y = y';
end

%Index must be between 2 and y_length-1
if (find(x==input) < 2 || find(x==input) > length(y)-1)
    error('Index is out of bounds. It must be between 2 and (signal_length)-1');
elseif isempty(find(x==input, 1))
    error('Enter an input that is on x');
else    
    index = find(x==input);
    surr_indices = [index-1;index;index+1];
    surr_x = x(surr_indices);
    A = [surr_x.^2 surr_x surr_x.^0];
    surr_y = y(surr_indices);
    coeffs = A\surr_y;
    x_loc = -coeffs(2)/(2*coeffs(1));
    interp_max = coeffs'*[x_loc.^2;x_loc;1];
end;
end