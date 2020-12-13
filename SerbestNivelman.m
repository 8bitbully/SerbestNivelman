% SerbestNivelman:
% Nivelman aglarinin serbest dengelenmesi,
% Odev icin yapılmistir.
%
% Dosyalari yuklemek icin SerbestNivelman.readFile() metodu kullanılmalı
% veya cell array şeklinde sınıfa verilmelidir.
%
% ölcü dosya yapısı:
% BN   SN       Yukseklik Farki(m)      Gecki Uzunlugu(km)
% {...}    {...}            {...}                               {...}
% {...}    {...}            {...}                               {...}
%   .        .               .                                   .
%   .        .               .                                   .
% yaklasik degerlerin dosya yapısı:
% Nokta No      Yaklasik Yukseklikler(m)
%     {...}                         {...}
%     {...}                         {...}
%       .                             .
%       .                             .
%

% version: 9.9.0.1467703 (R2020b)
% github.com/solounextracto
% @author: 
% @date: 20201204
classdef SerbestNivelman
    properties (Dependent = true, Access = public)
        A double
        l double
        P double
    end
    properties (Access = public)
        pointName(:, 1) cell
    end
    properties (Access = private)
        counter uint8
        measure
        elevation
        Measured
        f uint8
    end
    
    methods
        % constructor method
        function this = SerbestNivelman(measure, elevation, ct)
            if nargin < 3; this.counter = 1; else; this.counter = ct; end
            format longG
            this.measure = measure;
            this.elevation = elevation;
            this.Measured= this.measure;
        end
        
        function [EOA, params] = initAdjust(adjustmentObject)
            old_adjust = 0;
            while (adjustmentObject.counter - old_adjust) ~= 0
                [x, Qxx] = dengelemeBilinmeyen(adjustmentObject) ;
                [H, V] = kesinDeger(adjustmentObject, x) ;
                M = duyarlilik(adjustmentObject, V, Qxx) ;
                adjustmentObject = t_test(adjustmentObject, V, M) ;
            end
            EOA = adjustmentObject;
            params = struct('x', x, 'Qxx', Qxx, 'H', H, 'V', V, 'M', M);
        end
    end

    % SET METHODS
    methods
        % olculer
        function this = set.measure(this, msr)
            this.measure = msr ;
        end

        %yukseklikler
        function this = set.elevation(this, elev)
            this.elevation = elev ;
        end
        
        function this = set.Measured(this, ms)
            this.Measured =  ms;
        end
        
    end

    % GET METHODS 
    methods
        % serbestlik derecesi
        function f = get.f(this)
            n = size(this.Measured, 1) ;
            u = size(this.elevation, 1) ;
            d = 1 ; % Yükseklik ağlarının datum parametresi, 1 öteleme.
            
            f_ = n - u*1 + d ;
            f = f_ ;

            if f_ < 0
                error('varsayimlara dayali cozum yapilir.')
            elseif f_ == 0
                error('tek anlamli cebrik cozum yapilir.')
            end
        end
        
        function name = get.pointName(this)
            name = unique(cellfun(@num2str, this.Measured(:,1:2), 'UniformOutput', 0));
        end

        % stokastik model(agirlik matrisi)
        function p = get.P(this)
            p_ = [this.Measured{:, 4}];
            p_ = 1 ./ p_ ;
            p = diag(p_) ;
        end
         
        % katsayilar matrisi
        function a = get.A(this)
            A_ = zeros(size(this.Measured, 1), size(this.elevation, 1));
            points = this.pointName ;
            for i = 1 : size(this.Measured, 1)
               startidx = strcmp(string(points), string(this.Measured{i, 2}));
               stopidx = strcmp(string(points), string(this.Measured{i, 1}));
               A_(i, startidx) = 1;
               A_(i, stopidx) = -1;
            end
            a = A_;
        end
        
        % oteleme vektoru
        % TODO: try-catch-otherwise, l hesaplanamıyor ise exception fırlat.
        function l_ = get.l(this)
            elev = this.elevation ;
            msr = this.Measured ;
            try
                for i = 1 : size(msr, 1)
                    bn = elev{strcmp(string(elev(:,1)), string(msr{i, 2})), 2};
                    dn = elev{strcmp(string(elev(:,1)), string(msr{i, 1})), 2};
                    out.l(i) = bn - dn ;
                end
            catch
                error(' /// Lutfen nokta isimlerini kontrol ediniz.')
            end
            l__ = [msr{:, 3}] - out.l;
            l_ = round(l__'*1e3) ;
        end

    end
    
    methods
        % dengeleme bilinmeyeni hesabı
        %   "Tum iz minimum"
        %   args:
        %       this: object
        %   returns:
        %       x: serbest dengelenmis koordinatlar vektoru
        %       Qxx: Koordinatlari bilinmeyenlerin ters agirlik matrisi
        function [x, Qxx] = dengelemeBilinmeyen(this)
            A_ = this.A ;
            P_ = this.P ;
            l_ = this.l ;
            
            N = A_' * P_ * A_ ;
            n = A_' * P_ * l_ ;
            
            % pinv()
            G = 1./sqrt(size(A_, 2)) * ones(size(A_, 2), 1) ;
            E = G*G' ;
            
            Np = (N + E)^-1 - E ; % pseudo ters
            Qxx_ = Np ;
            Qxx = Qxx_ ;

            x = Qxx_ * n ;
        end

        % kesin değer hesabı
        %   args:
        %       this: Object
        %       x: serbest dengelenmis koordinatlar vektoru
        %   returns:
        %       H: Nokta yuksekliklerinin kesin degeri
        %       V: Duzeltmeler
        function [H, V] = kesinDeger(this, x)
            h = [this.elevation{:, 2}]' ;
            A_ = this.A ;
            l_ = this.l ;
            H = h + x / 1e3 ;
            V = A_*x - l_ ;
        end

        % duyarlılık hesapları
        %   args:
        %       V: Duzeltmeler
        %       Qxx: Bilinmeyenlerin ters agirlik matrisi
        %   returns:
        %       M: Duyarlilik hesaplari
        function M = duyarlilik(this, V, Qxx)
            A_ = this.A ;
            P_ = this.P ;
            f_ = this.f ;
            Qxx_ = Qxx ;

            m0 = sqrt((V'*P_*V) / (f_)) ; % birim olcunun ortalama hatasi
            mi = m0 ./ sqrt(diag(P_)) ; % olculerin ortalama hatasi
            mx = m0 .* sqrt(diag(Qxx_)) ; % bilinmeyenlerin ortalama hatasi
            Qll_ = A_ * Qxx_ * A_' ; % dengeli olculerin ters agirlik matrisi
            ml_ = m0 * sqrt(diag(Qll_)) ; % dengeli olculerin ortalama hatasi
            Qvv = P_^-1 - Qll_ ; % duzeltmelerin tes agirlik matrisi
            mv = m0 * sqrt(diag(Qvv)) ; % duzeltmelerin ortalama hatasi

            M = struct( ...
                'm0', m0, ...
                'mi', mi, ...
                'mx', mx, ...
                'Qll', Qll_, ...
                'ml', ml_, ...
                'Qvv', Qvv, ...
                'mv', mv) ;
        end

        % Uyusumsuz olculer testi, 
        %   args:
        %       V: Duzeltmeler
        %       M: Duyarlilik hesaplari
        %       alpha: 
        %   returns:
        %       Object: SerbestDengeleme nesnesi
        function Object = t_test(this, V, M, alpha)
            if nargin < 4; alpha = 5; end
            Qvv = M.Qvv ;
            P_ = this.P ;
            f_ = this.f ;

            n = this.counter;

            s0 = sqrt((1 / (f_ - 1)) * ((V'*P_ * V) - (V.^2) ./ (diag(Qvv)))) ;

            T = abs(V) ./ (s0.*sqrt(diag(Qvv))) ;

            T_ = max(T) ;
            t_ = tinv(1 - (alpha/1e2)/2, f_ - 1) ;

            fprintf('// %d. kez dengeleme yapilmistir.\n', n) ;

            if T_ > t_
                fprintf('// "T:%f > t:%f" oldugundan\n', T_, t_)
                fprintf('// %d. olcu uyusumsuzdur ve olculerin arasindan cikarilmistir.\n', find(T==T_)) ;
                idx = ( T == T_ ) ;
                this.Measured(idx, :) = [];
                n = n + 1;
            else
                fprintf('// "T:%f < t:%f" oldugundan\n', T_, t_)
                fprintf('// uyusumsuz olcu yoktur.\n')
                n = 0;
            end
            fprintf('%s\n', char(ones(1,100)*45))
            Object = SerbestNivelman(this.Measured, this.elevation, n) ;
        end
    end
    
    methods (Static)
        % Dosya Okuma: only .txt
        %   args:
        %       filename: str
        %   returns:
        %       content: cell array
        function content = readFile(filename)
            data = readcell(filename) ;
            check = cellfun(@ismissing, data, 'UniformOutput', false);
            ret = cellfun(@all, check);
            count = 1; err.c = [];
            for i = 1:numel(ret)
                if ret(i); err.c(count) = mod(i, size(ret, 1)); count = count + 1; end
            end
            err = unique(err.c);
            if numel(err) >= 1
                error([mat2str(err), ': satirlari hatalidir.'])
            end
            content = data;
        end
    end
end