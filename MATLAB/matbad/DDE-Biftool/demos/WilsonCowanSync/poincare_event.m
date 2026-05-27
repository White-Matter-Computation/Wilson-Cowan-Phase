function [value, isterminal, direction] = poincare_event(~, y, ~, u_star)
    value = y(1) - u_star;

    % Do not stop integration at events.
    isterminal = 0;

    % Detect upward crossings only.
    direction = 1;
end