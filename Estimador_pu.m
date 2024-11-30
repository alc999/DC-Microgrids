classdef Estimador_pu < matlab.System

    % Public, tunable properties
    properties
        Vbase = 380;
        Pbase = 1000;
    end

    properties (DiscreteState)

    end

    % Pre-computed constants
    properties (Access = private)
        V = zeros(4, 20)
        I = zeros(4, 20)
        paso = 5e-5;
        ni = 10;
%         paso = 0.025;
        gL = zeros(4,4);
        gcc = zeros(4,4);
        goo = zeros(4,4);
        gco = zeros(4,4);
        Q = (1000/380)*eye(4,4);
        Gk = zeros(4,4)
        dv = 1e-6;
        Rbase;
    end

    methods (Access = protected)
        function setupImpl(obj)
            Ibase = obj.Pbase/obj.Vbase;
            obj.Rbase = obj.Vbase/Ibase;
            %R = [0.021, 0.0525, 0.0105, 0.11462, 0.0168, 0.042, 0.08336];
%             R = [0.0766, 0.087, 0.1704, 0.1044, 0.1136, 0.0557, 0.1392];
            r = [0.0226 0.0257 0.0769 0.0164 0.0513 0.0164 0.0411];
            R = r/obj.Rbase;
            obj.gcc = [((1/R(1))+(1/R(3))),     0,      0,       0;
                0,           1/R(2),    0,       0;
                0,              0,    1/R(4),    0;
                0,              0,      0,       1/R(6)];
                    
            obj.gco = [-1/R(1), -1/R(3),    0,      0;
                   -1/R(2),    0,       0,      0;
                      0,    -1/R(4),    0,      0;
                      0,       0,    -1/R(6),   0];
                  
                  
            obj.goo = [(1/R(1))+(1/R(2)),              0,                          0,                      0;
                               0,               (1/R(3))+(1/R(4))+(1/R(5)),          -1/R(5),                   0;
                               0,                       -1/R(5),            (1/R(5))+(1/R(6))+(1/R(7)),      -1/R(7);
                               0,                           0,                       -1/R(7),           (1/R(7))];

            obj.Gk = obj.gcc - obj.gco*((obj.goo+obj.gL)\obj.gco');
        end

        function Gbus = stepImpl(obj,vk, pk)
            % Implement algorithm. Calculate y as a function of input u and
            % discrete states.
            obj.V = circshift(obj.V, [0 1]);
            obj.V(:,1) = vk/obj.Vbase;
            obj.I = circshift(obj.I, [0 1]);
            obj.I(:,1) = (pk/obj.Pbase)./(vk/obj.Vbase);

            
            y = obj.gcc*obj.V-obj.I;
%             for iter = 1:50
%                 t1 = (obj.gL+obj.goo)\(obj.gco'*obj.V);
%                 gradiente = (obj.gL+obj.goo)'\(obj.gco'*(y-obj.gco*t1)*t1');
%                 obj.gL = obj.gL - obj.paso*gradiente;
%             end

            % Enfoque con cuadratica
            if (obj.V(1,1) > (obj.V(1,2) - obj.dv)) && (obj.V(1,1) < (obj.V(1,2) + obj.dv))
                Gbus = (obj.Gk+obj.Gk')/2;
            else
                for iter = 1:obj.ni
                    z = (obj.gL+obj.goo)\(obj.gco'*obj.V);
                    gradiente = (obj.gL+obj.goo')\(obj.gco'*obj.Q*(y-obj.gco*z)*z');
                    obj.gL = obj.gL - obj.paso*gradiente;
                end
                obj.Gk = obj.gcc - obj.gco*((obj.goo+obj.gL)\obj.gco');
                Gbus = ((obj.Gk+obj.Gk')/2);
            end
         
        end

        function resetImpl(obj)
            % Initialize / reset discrete-state properties
        end
    end
end
