clearvars; clc; 

serbest = SerbestNivelman('data/olcu.txt', 'data/yaklasik.txt');

[x1, Qxx] = dengelemeBilinmeyen(serbest) ;
[H1, V] = kesinDeger(serbest, x1) ;
M1 = duyarlilik(serbest, V, Qxx) ;
obj = t_test(serbest, V, M1.Qvv) ;

[x2, Qxx] = dengelemeBilinmeyen(obj) ;
[H2, V] = kesinDeger(obj, x2) ;
M2 = duyarlilik(obj, V, Qxx) ;
obj2 = t_test(obj, V, M2.Qvv) ;
