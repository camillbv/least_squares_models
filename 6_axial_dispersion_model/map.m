function [x_desired,w_desired] = map(x_current,w_current,des_min,des_max)
%% header
% author(s)             :: c.b.vik
% date                  :: 27.06.2013 (revised)
% purpose               :: map a set of points and weights into desired 
%                          domain from current domain
%
%% pack out domain boundaries from input
x_cur_min = x_current(1);
x_cur_max = x_current(length(x_current));
x_des_min = des_min;
x_des_max = des_max;
%
%% calculate domain lengths 
delta_cur = x_cur_max - x_cur_min; % size of current domain
delta_des = x_des_max - x_des_min; % size of desired domain
%
%% calculate determinant
det = delta_des/delta_cur;
%
%% main
x_desired = det*(x_current-x_cur_min)+x_des_min; % map points
w_desired = det*w_current;                       % map weights
% 
