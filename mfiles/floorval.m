function floorval = floorval(val,nS)
    pw = ceil(log10(abs(val))); % Find order of magnitude of val.
    res = 10^(pw-nS); % Resolution to round to.
    floorval = floor(val/res)*res; % < change floor() to ceil(), for ceiling equivalent.
end
