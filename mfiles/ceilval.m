function ceilval = ceilval(val,nS)
    pw = ceil(log10(abs(val))); % Find order of magnitude of val.
    res = 10^(pw-nS); % Resolution to round to.
    ceilval = ceil(val/res)*res; % < change ceil() to floor(), for flooring equivalent.
end
