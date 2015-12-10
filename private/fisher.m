function z = fisher(r)
% Fisher's Z-transform.
% z = fisher(r) returns the Fisher's Z-transform of the correlation
% coefficient r

z = 0.5 * log((1 + r) ./ (1 - r));