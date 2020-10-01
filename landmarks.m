function f = landmarks()
    step = 0.4;
    n1 = 4/step;
    n2 = n1*2;
    n3 = n1*3;
    f = zeros(2,n3);
    for i = 1:n1
        f(1,i) = -2+step*(i);
        f(2,i) = -1;
    end
    for i = 1:n1
        f(1,i+n1) = -2+0.5*step*(i-1);
        f(2,i+n1) = (sqrt(3)+1)/2 * (-2+0.5*step*(i-1)) + sqrt(3);
    end
    for i = 1:n1
        f(1,i+n2) = 0.5*step*(i-1);
        f(2,i+n2) = (-(sqrt(3)+1)/2 * (0.5*step*(i-1)) + sqrt(3));
    end
end
