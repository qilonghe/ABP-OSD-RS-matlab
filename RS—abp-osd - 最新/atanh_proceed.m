function y = atanh_proceed(x)
y =atanh(x);
if y == Inf
    y = 100000;
end
if y == -Inf
    y = -100000;
end
end