function [W,R] = Unpack(S)
W = S(1:6);
R = [S(7:9),S(10:12),S(13:15)];
end