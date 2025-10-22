%Utility function:

function u_func = u(c, gamma)
    if gamma == 1
       % u_func = log(c);
    else
        u_func = (c.^(1 - gamma)) ./ (1 - gamma);
    end
end
