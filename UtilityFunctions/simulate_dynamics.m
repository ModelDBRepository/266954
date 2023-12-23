function [r, v] = simulate_dynamics(I, N, T, w, dt, tau)
    tau_s = 5;
    r = zeros(N, length(T));
    v = zeros(N, length(T));
    for i = 1:length(T)-1
        inp = w' * r(:,i) + I(:,i);
        v(:,i+1) = v(:,i) + dt/tau_s * (-v(:,i)+inp);
        rd = -r(:,i) + NL(inp);
        r(:,i+1) = r(:,i) + dt/tau * rd;
    end
end

function [zo] = NL(zi)
    zo = zi;
    zo(zo < 0) = 0;
    zo = zo.^1;
end