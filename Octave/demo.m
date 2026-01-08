%Only needed once:
%pkg install -forge symbolic-3.0.1.tar.gz
pkg load symbolic

syms x;
e0 = sym([x-1, x, 0]);
n=3;
d=1;

Invariants = invariants(n, d, @substitutionSumtwothree, e0);
