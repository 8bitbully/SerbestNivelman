% version: 9.9.0.1467703 (R2020b)
%
clc;clear;clearvars;
mainfun = localfunctions;


% // dosyadan olculeri okuma
olcu = SerbestNivelman.readFile('foy-olcu.txt');
yaklas = SerbestNivelman.readFile('foy-yaklasik.txt');

% // dengeleme icin nesnenin olusturulmasi
serbest = SerbestNivelman(olcu, yaklas);
% // en son dengeleme degerlerini verir.
[EOA, params] = initAdjust(serbest);

% 1. dengeleme
% [x1, Qxx1] = dengelemeBilinmeyen(serbest) ;
% [H1, V1] = kesinDeger(serbest, x1) ;
% M1 = duyarlilik(serbest, V1, Qxx1) ;
% obj1 = t_test(serbest, V1, M1) ;