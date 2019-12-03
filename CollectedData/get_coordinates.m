function coord = get_coordinates(location, x_v, y_v, z_v)

% This function looks complicated, but it's just the same as
% get_coordinates_old.m, except that I add the extreme case when the x, y
% or z coordinate could be the length()

dim_z = floor(location/((length(y_v))*length(x_v)));
zresd = location - dim_z*length(y_v)*length(x_v);
if zresd == 0
    z_coordinate = z_v(dim_z);
    y_coordinate = y_v(length(y_v));
    x_coordinate = x_v(length(x_v));
    coord = [x_coordinate, y_coordinate, z_coordinate];
else
    z_coordinate = z_v(dim_z + 1);
    dim_y = floor((location - dim_z*length(y_v)*length(x_v))/length(x_v));
    yresd = location - dim_z*length(y_v)*length(x_v) - dim_y*length(x_v);
    if yresd == 0
        y_coordinate = y_v(dim_y);
        x_coordinate = x_v(length(x_v));
        coord = [x_coordinate, y_coordinate, z_coordinate];
    else
        y_coordinate = y_v(dim_y + 1);
        x_coordinate = x_v(yresd);
        coord = [x_coordinate, y_coordinate, z_coordinate];
    end
end

end
