function J = JacobianEst(func, x0, opts)

arguments
    func % function handle func(x)
    x0 % location to estimate Jacobian
    opts.step
    opts.epsilon = 1e-3;
    opts.xlb = []
    opts.xub = []

end


% Number of inputs
Nx = length(x0);

% Evaluate function at x0
y0 = func(x0);

% Number of outputs
Ny = length(y0);

% Initialise Jacobian matrix
J = zeros(Ny, Nx);


% Find partial derivative for each input element
for xindx = 1:Nx

    dx = zeros(size(x0));
    thisstep = opts.step(xindx);

    % Adjust step to respect bounds
    if ~isempty(opts.xub) && x0(xindx) + thisstep > opts.xub(xindx)
        thisstep = -thisstep;
    elseif ~isempty(opts.xlb) && x0(xindx) - thisstep < opts.xlb(xindx)
        thisstep = -thisstep;
    end


    dx(xindx) = thisstep;
    J(:, xindx) = (func(x0 + dx/2) - func(x0-dx/2)) / thisstep;

end





end