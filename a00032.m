clc;clear;clearvars;
main = localfunctions;
% // a00032.m scripti run(F5) ile çalıştırıldığında sonuçları
% // verir ve excele yazdırır

% // eger ayni dosya varsa siler.
prefix = 'Denge-';
for i = 1:3
    if isfile([prefix,num2str(i),'.xlsx'])
        delete([prefix,num2str(i),'.xlsx'])
    end
end


% // dosyadan olculeri okuma
olcu = SerbestNivelman.readFile('olculer.txt');
yaklas = SerbestNivelman.readFile('yaklasikyukseklik.txt');

% // dengeleme icin nesnenin olusturulmasi
serbest = SerbestNivelman(olcu, yaklas);
% // en son dengeleme degerlerini verir.
% [EOA, params] = initAdjust(serbest);

% 1. dengeleme
[x1, Qxx1] = dengelemeBilinmeyen(serbest) ;
[H1, V1] = kesinDeger(serbest, x1) ;
M1 = duyarlilik(serbest, V1, Qxx1) ;
obj1 = t_test(serbest, V1, M1) ;
m0 = M1.m0; mx = M1.mx; l = serbest.l; A = serbest.A; Qvv = M1.Qvv;
% // sonuclarin excel e yazdirilmasi
writeExcel('Denge-1', x1, Qxx1, H1, V1, m0, mx, A, l, Qvv)


% // Functions
function writeExcel(filename, varargin)
    n = length(varargin);
    filename = [filename, '.xlsx'];
    for i = 1 : n
        writematrix(varargin{i}, filename,'Sheet',inputname(i+1));
    end
end