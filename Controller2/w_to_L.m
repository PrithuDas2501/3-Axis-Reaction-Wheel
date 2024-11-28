function L = w_to_L(w)

L1 = diag(w);

L2 = [0, w(3), w(2);
      w(3), 0, w(1);
      w(2), w(1),0];

L = [L1,L2];
end