function y = plotneuralnetwork(weights)

%u = -10.00:0.01:10.00;
u = 1.4:0.001:3.6;

reg = neuralnetwork(weights, u);

plot(u,reg);

end