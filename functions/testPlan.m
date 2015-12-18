clear;
load('control.mat')
ec = ec';
ed = ed';

[ee_test, eff]  = plan3gs(ec, ed);
ef_test = plan3gf(ec,t,ed,es);
Ke_test = plan3ge(ec,t,D,ed,es);

Ke_diff_norm = norm(Ke-Ke_test)
ef_diff_norm = norm(ef'-ef_test)
es_diff_norm = norm(D*ee_test - es)