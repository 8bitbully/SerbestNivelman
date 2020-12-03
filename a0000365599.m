clc;clear;clearvars;
main = localfunctions;
prefix = 'Dengeleme-';
for i = 1:3
    if isfile([prefix,num2str(i),'.xlsx'])
        delete([prefix,num2str(i),'.xlsx'])
    end
end


% dosyadan olculeri okuma
olcu = SerbestNivelman.readFile('olcu.txt');
yaklas = SerbestNivelman.readFile('yaklasik.txt');

% dengeleme icin nesnenin olusturulmasi
serbest = SerbestNivelman(olcu, yaklas);

% 1. dengeleme
[x1, Qxx1] = dengelemeBilinmeyen(serbest) ;
[H1, V1] = kesinDeger(serbest, x1) ;
M1 = duyarlilik(serbest, V1, Qxx1) ;
obj1 = t_test(serbest, V1, M1) ;
m0 = M1.m0; mx = M1.mx; l = obj1.l;
% sonuclarin excel e yazdirilmasi
writeExcel('Dengeleme-1', x1, Qxx1, H1, V1, m0, mx, l)

% 2. dengeleme
[x2, Qxx2] = dengelemeBilinmeyen(obj1) ;
[H2, V2] = kesinDeger(obj1, x2) ;
M2 = duyarlilik(obj1, V2, Qxx2) ;
obj2 = t_test(obj1, V2, M2) ;
m0 = M1.m0; mx = M1.mx; l = obj1.l;
% sonuclarin excel e yazdirilmasi
writeExcel('Dengeleme-2', x2, Qxx2, H2, V2, m0, mx, l)

% 3. dengeleme
[x3, Qxx3] = dengelemeBilinmeyen(obj2) ;
[H3, V3] = kesinDeger(obj2, x3) ;
M3 = duyarlilik(obj2, V3, Qxx3) ;
obj3 = t_test(obj2, V3, M3) ;
m0 = M1.m0; mx = M1.mx; l = obj1.l;
% sonuclarin excel e yazdirilmasi
writeExcel('Dengeleme-3', x3, Qxx3, H3, V3, m0, mx, l)

% Functions
function writeExcel(filename, varargin)
    n = length(varargin);
    filename = [filename, '.xlsx'];
    for i = 1 : n
        writematrix(varargin{i}, filename,'Sheet',inputname(i+1));
    end
end