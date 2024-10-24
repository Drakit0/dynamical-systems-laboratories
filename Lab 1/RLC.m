function [ft1, ft2] = RLC(R, L, C, mode)
% Calculate:
% - ft_1: V0(s)/Vi(s)
% - ft2: I(s)/Vi(s)
% of the given circuit through 3 diff methods 
    
    if mode == "tf" % Polynomial division
        ft1 = tf([1],[C*L, C*R, 1]);
        ft2 = tf([C,0],[C*L, C*R, 1]);

    else if mode == "zpk" % Zeros, poles and static gain
        z_1 = [];
        p_1 = -R/(2*L);
        k_1 = 1/(C*L);

        z_2 = 0;
        p_2 = p_1;
        k_2 = 1/L;

        ft1 = zpk(z_1, [p_1, p_1], k_1);
        ft2 = zpk(z_2, [p_2, p_2], k_2);

    else % Symbolic s operations 
        s = tf('s');

        ft1 = 1/(C*L*s^2 + C*R*s + 1);
        ft2 = C*s/(C*L*s^2 + C*R*s + 1);

    end

end

