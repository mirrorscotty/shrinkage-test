classdef data
    %DATA Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %Bulk modulus of material
        k1 = 500;
        k2 = 400;
        muk = 3;
        
        %Bulk modulus of solid only
        ks1 = 600;
        ks2 = 500;
        muks = 2;
        
        %Shear modulus
        g1 = 500;
        g2 = 350;
        mug = 1;
        
        %Bulk modulus of fluid
        kf = 300;
        
        %Porosity
        poro = .6;
        
        %Densities of both solid and fluid
        rhos = 2000;
        rhof = 1000;
        
        %Permutivity
        kappa = 42;
    end
    
    methods
        function [ output_args ] = alphahat(this, s)
            output_args = 1 - this.khat(s)/this.kshat(s);
        end

        function [ output_args ] = beta(this,s)
            output_args = this.poro^2*s*this.kappa*this.rhof/ ...
            (this.poro^2+s*this.kappa(this.rhoa(s)+this.poro*this.rhof));
        end

        function [ output_args ] = ghat(this,s)
            g = this.g1*this.g2/(this.g1+this.g2);
            qg = this.mug/this.g2;
            pg = this.mug/(this.g1+this.g2);

            output_args = g*(1+qg*s)/(1+pg*s);
        end

        function [ output_args ] = khat( this,s )
            k = this.k1*this.k2/(this.k1+this.k2);
            qk = this.muk/this.k2;
            pk = this.muk/(this.k1+this.k2);

            output_args = k*(1+qk*s)/(1+pk*s);
        end

        function [ output_args ] = kshat( this,s )
            ks = this.ks1*this.ks2/(this.ks1+this.ks2);
            qks = this.muks/this.ks2;
            pks = this.muks/(this.ks1+this.ks2);

            output_args = ks*(1+qks*s)/(1+pks*s);
        end

        function [ output_args ] = rhat( this,s )
            ks = this.kshat(s);
            kff = this.kf; %Renamed because Matlab likes to warn about pointless things.
            k = this.khat(s);
            p = this.poro;

            output_args = p^2*kff*ks/(kff*(ks-k)+p*ks*(ks-kff));
        end

        function [ output_args ] = rho(this,s)
             output_args = this.rhos*(1-this.poro) + this.poro*this.rhof;
        end

        function [ output_args ] = rhoa( this,s )
            output_args = 0.66*this.poro*this.rhof;
        end
        
    end
    
end

