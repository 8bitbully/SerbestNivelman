classdef SerbestNivelman
    properties
        A double
        l double
        P double
    end
    properties (Access = private)
        measure
        elevation
        Measured
        f double
        pointName(:, 1) cell
    end

    % constructor method
    methods
        function this = SerbestNivelman(measure, elevation)
            this.measure = measure;
            this.elevation = elevation;
            this.Measured = this.measure ;
        end
    end

    % SET METHODS
    methods
        function this = set.measure(this, msr)
            this.measure = readcell(msr) ;
        end
        
        function this = set.elevation(this, elev)
            this.elevation = readcell(elev) ;
        end
        
        function this = set.Measured(this, ms)
            this.Measured =  ms;
        end
        
    end

    % GET METHODS 
    methods
        function f = get.f(this)
            n = size(this.Measured, 1) ;
            u = size(this.elevation, 1) ;
            d = 1 ;
            
            f_ = n - u*1 + d ;
            f = f_ ;

            if f_ < 0
                error('dengeleme yapilmaz')
            elseif f_ == 0
                error('dengeleme == ?')
            end
        end

        function p = get.P(this)
            p_ = [this.Measured{:, 4}];
            p_ = 1 ./ p_ ;
            p = diag(p_) ;
        end
        
        function a = get.A(this)
            A_ = zeros(size(this.Measured, 1), size(this.elevation, 1));
            points = unique(cellfun(@num2str, this.Measured(:,1:2), 'UniformOutput', 0));
            for i = 1 : size(this.Measured, 1)
               startidx = strcmp(string(points), string(this.Measured{i, 2}));
               stopidx = strcmp(string(points), string(this.Measured{i, 1}));
               A_(i, startidx) = 1;
               A_(i, stopidx) = -1;
            end
            a = A_;
        end
        
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
                error(' /// Lütfen nokta isimlerini kontrol ediniz.')
            end
            l__ = [msr{:, 3}] - out.l;
            l_ = round(l__'*1e3) ;
        end

    end
    
    methods
        % dengeleme bilinmeyeni hesabı
        function [x, Qxx] = dengelemeBilinmeyen(this)
            A_ = this.A ;
            P_ = this.P ;
            l_ = this.l ;
            
            N = A_' * P_ * A_ ;
            n = A_' * P_ * l_ ;
            
            G = 1./sqrt(size(A_, 2)) * ones(size(A_, 2), 1) ;
            E = G*G' ;
            
            Np = (N + E)^-1 - E ;
            Qxx_ = Np ;
            Qxx = Qxx_ ;

            x = Qxx_ * n ;
        end

        % kesin değer hesabı
        function [H, V] = kesinDeger(this, x)
            h = [this.elevation{:, 2}]' ;
            A_ = this.A ;
            l_ = this.l ;
            H = h + x / 1e3 ;
            V = A_*x - l_ ;
        end

        % duyarlılık hesapları
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
            Qvv = P_^-1 - Qll_ ;
            mv = m0 * sqrt(diag(Qvv)) ;

            M = struct( ...
                'm0', m0, ...
                'mi', mi, ...
                'mx', mx, ...
                'Qll', Qll_, ...
                'ml', ml_, ...
                'Qvv', Qvv, ...
                'mv', mv) ;
        end

        % Uyuşumsuz ölçüler testi, 
        % <output>Object: SerbestNivelman</output>
        function testObject = t_test(this, V, Qvv, alpha)
            if nargin < 4; alpha = 5; end
            P_ = this.P ;
            f_ = this.f ;

            persistent n
            if isempty(n); n = 0; end
            
            n = n + 1 ;

            s0 = sqrt((1 / (f_ - 1)) * ((V'*P_ * V) - (V.^2) ./ (diag(Qvv)))) ;

            T = abs(V) ./ (s0.*sqrt(diag(Qvv))) ;

            T_ = max(T) ;
            t_ = tinv(1 - (alpha/1e2)/2, f_ - 1) ;

            fprintf('->%d. kez dengeleme yapilmistir.\n', n) ;

            if T_ > t_
                fprintf('->%d. olcu uyusumsuzdur ve olcularin arasindan cikarilmistir.\n', find(T==T_)) ;
                idx = ( T == T_ ) ;
                this.Measured(idx, :) = [];
            else
                fprintf('->uyusumsuz olcu yoktur.\n')
                n = 0;
            end
            fprintf('----------------------\n')
            testObject = this ;
        end
    end
end