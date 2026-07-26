function y=my_round(x)
if x<0
    y=-ceil(abs(x));
else
    y=ceil(x);
end