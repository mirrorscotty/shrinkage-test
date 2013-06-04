classdef poroelastic_analytic
    % Analytic solution to a 1D poroelastic column under stress
    
    properties
        p0 = 1;
        l = 1;
        k = 3e7;
        ks = 3.6e10;
        kf = 2.3e9;
        g = 12.2e7;
        phi = .79;
        kappa = 1e-8;
        rhof = 1000;
        rhos = 1396;
    end
    
    methods
        function [out] = Ehat(this, s)
            out = this.k + 4/3*this.g;
        end
        function [out] = rhoa(this)
            out = 0.66*this.phi*this.rhof;
        end
        function [out] = rho(this)
            out = this.rhos*(1-this.phi)+this.phi*this.rhof;
        end
        function [out] = alpha(this, s)
            out = 1+this.k/this.ks;
        end
        function [out] = beta(this, s)
            p = this.phi;
            K = this.kappa;
            ra = this.rhoa;
            rf = this.rhof;
            out = p^2*s*K*rf/(p^2+s*K*(ra+p*rf));
        end
        function [out] = Rhat(this, s)
            p = this.phi;
            kf = this.kf;
            ks = this.ks;
            k = this.k;
            out = p^2*kf*ks^2/(kf*(ks-k)+p*ks*(ks-kf));
        end
        function [out] = Aq(this, s)
            out = this.Ehat(s)*this.phi^2/this.rhof;
        end
        function [out] = Bq(this, s)
            E = this.Ehat(s);
            p = this.phi;
            R = this.Rhat(s);
            A = this.alpha(s);
            B = this.beta(s);
            r = this.rho();
            rf = this.rhof();
            out = E*p^2/R + (r-B*rf)*B/rf+(A-B)^2;
        end
        function [out] = Cq(this, s)
            out = this.phi^2*(this.rho() - this.beta(s)*this.rhof)/this.Rhat(s);
        end
        function [out] = root1(this, s)
            Aq = this.Aq(s);
            Bq = this.Bq(s);
            Cq = this.Cq(s);
            out = sqrt((Bq+sqrt(Bq^2-4*Aq*Cq))/(2*Aq));
        end
        function [out] = root2(this, s)
            Aq = this.Aq(s);
            Bq = this.Bq(s);
            Cq = this.Cq(s);
            out = sqrt((Bq-sqrt(Bq^2-4*Aq*Cq))/(2*Aq));
        end
        function [out] = di(this, root, s)
            if(root == 1)
                r = this.root1(s);
            else
                r = this.root2(s);
            end
            E = this.Ehat(s);
            rho = this.rho();
            rhof = this.rhof;
            A = this.alpha(s);
            B = this.beta(s);
            out = (E*r^2-(rho-B*rhof))/((A-B)*r);
        end
        function [out] = uhaty(this, s, y)
            E = this.Ehat(s);
            d1 = this.di(1, s);
            d2 = this.di(2, s);
            r1 = this.root1(s);
            r2 = this.root2(s);
            p = this.p0;
            l = this.l;
            
            out = p/(E*(d1*r2-d2*r1))* ...
                (d2*(exp(-r1*s*(l-y))-exp(-r1*s*(l+y))) / (s*(1+exp(-2*r1*s*l)))- ...
                 d1*(exp(-r2*s*(l-y))-exp(-r2*s*(l+y))) / (s*(1+exp(-2*r2*s*l))));
        end
        
        function [out] = plotu(this, start, fin)
            s = linspace(start, fin);
            for i = 1:length(s)
                u(i) = this.uhaty(s(i), 1);
            end
            out = u;
            plot(s,u);
        end
    end
    
end

