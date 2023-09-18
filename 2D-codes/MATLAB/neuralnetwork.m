function y = neuralnetwork(weights,u)

L = length(weights);
res = u;

for i=1:(L-1)

    res = (4/3.14159267)*atan(weights(i)*res);%tanh(weights(i)*res);
    
end

y = weights(L)*res;

end