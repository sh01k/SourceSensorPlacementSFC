function [pos] = rect_perim(num,len_x,len_y,shift_x,shift_y)
% Discretize perimeter of rectangular region at regular intervals
% INPUT
%    num              Number of positions
%    len_x, len_y     Length in x and y
%    shift_x, shift_y Shift length in x and y from origin
% OUTPUT
%    pos              Position vectors
%
% Jun 2019 Shoichi Koyama, Gilles Chardon, and Laurent Daudet

%Interval
d = (len_x*2+len_y*2)/num;
perim_vec = (((1:num)-0.5)*d).';

perim_h_u = perim_vec(perim_vec<=len_x);
x = perim_h_u-len_x/2+shift_x;
y = (len_y/2+shift_y)*ones(length(perim_h_u),1);

perim_v_r = perim_vec(perim_vec>len_x & perim_vec<=len_x+len_y)-len_x;
x = [x; (len_x/2+shift_x)*ones(length(perim_v_r),1)];
y = [y; -perim_v_r+len_y/2+shift_y];

perim_h_d = perim_vec(perim_vec>len_x+len_y & perim_vec<=len_x*2+len_y)-len_x-len_y;
x = [x; -perim_h_d+len_x/2+shift_x];
y = [y; (-len_y/2+shift_y)*ones(length(perim_h_d),1)];

perim_v_l = perim_vec(perim_vec>len_x*2+len_y & perim_vec<=len_x*2+len_y*2)-len_x*2-len_y;
x = [x; (-len_x/2+shift_x)*ones(length(perim_v_l),1)];
y = [y; perim_v_l-len_y/2+shift_y];

pos = [x, y];

end