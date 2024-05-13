function [BL_MESH] = Set_Grid_new(BL_MESH)
N_BL=BL_MESH.N_BL;
ymax_BL=BL_MESH.ymax_BL;
yi_BL=BL_MESH.yi_BL;

% Chebyshev Configuration Point
z = -cos((0:(N_BL-1))/(N_BL-1)*pi);

a = ymax_BL * yi_BL / (ymax_BL - 2.0 * yi_BL);
b = 1 + 2.0 * a / ymax_BL;  

% Coordinate Transformation
y_BL = zeros(N_BL,1); dzdy_BL = zeros(N_BL,1); d2zdy2_LST = zeros(N_BL,1);

    for i=1:N_BL
        y_BL(i) = a * (1 + z(i)) / (b - z(i));
        dzdy_BL(i) = a * (1 + b) / (a + y_BL(i))^2;
        d2zdy2_BL(i) = - 2 * a * (1 + b) / (a + y_BL(i))^3;
    end
    
BL_MESH.y_BL=y_BL;
BL_MESH.dzdy_BL=dzdy_BL;
BL_MESH.d2zdy2_LST=d2zdy2_BL;

[D_cheb] = chebDiff(N_BL);

Dy_BL = diag(dzdy_BL) * D_cheb; %1st order
DDy_BL= Dy_BL* Dy_BL; %2nd order

BL_MESH.Dy_BL=Dy_BL;
BL_MESH.DDy_BL=DDy_BL;
end

