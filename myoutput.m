function stop = myoutput(x, optimValues, state)
stop = false;
%global x_opt;
%x_opt(optimValues.iteration+1,:)=x;
end