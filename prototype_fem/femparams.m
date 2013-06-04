classdef femparams
    properties
        %Time stuff
        tinit = 0;
        tend = 100;
        dt = 1;
        
        %Space stuff
        xinit = 0;
        xend = 1;
        dx = .5;
        
        %Initial condition
        ICu = 0;
        ICp = 1;
    end
    
    methods
        function [z] = phi0(x)
            z = 1-x;
        end
        function [z] = phi1(x)
            z = x;
        end
        function [z] = dphi0(x)
            z = -1;
        end
        function [z] = dphi1(x)
            z = 1;
        end
    end
end

