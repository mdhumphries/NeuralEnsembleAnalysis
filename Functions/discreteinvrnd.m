function X = discreteinvrnd(p,m,n)

% p = vector of discrete probabilities; m,n = size of returned matrix

X = zeros(m,n); % Preallocate memory
cum_p = cumsum(p);
for i = 1:m*n
    u = rand;
    I = find(u < cum_p);
    try
    X(i) = min(I);
    catch
        keyboard
    end
end
