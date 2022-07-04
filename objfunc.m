%Name - Mridul Gupta
%Roll No - 19IM30025

function cost = objfunc(x,y,rD,D)
    cost = sum(x.*y)-rD*D;
    %cost = sum(x.*y);
end                                               