function [LST_MESH] = Set_LST_Grid_new(LST_MESH)
N_LST=LST_MESH.N_LST;
ymax_LST=LST_MESH.ymax_LST;
yi_LST=LST_MESH.yi_LST;

% Chebyshev Configuration Point
z = -cos((0:(N_LST-1))/(N_LST-1)*pi);

a = ymax_LST * yi_LST / (ymax_LST - 2.0 * yi_LST);
b = 1 + 2.0 * a / ymax_LST;  

% Coordinate Transformation
y_LST = zeros(N_LST,1); dzdy_LST = zeros(N_LST,1); d2zdy2_LST = zeros(N_LST,1);

    for i=1:N_LST
        y_LST(i) = a * (1 + z(i)) / (b - z(i));
        dzdy_LST(i) = a * (1 + b) / (a + y_LST(i))^2;
        d2zdy2_LST(i) = - 2 * a * (1 + b) / (a + y_LST(i))^3;
    end
    
LST_MESH.y_LST=y_LST;
LST_MESH.dzdy_LST=dzdy_LST;
LST_MESH.d2zdy2_LST=d2zdy2_LST;

[D_cheb] = chebDiff(N_LST);

Dy_LST = diag(dzdy_LST) * D_cheb; %1st order
DDy_LST= Dy_LST* Dy_LST; %2nd order

LST_MESH.Dy_LST=Dy_LST;
LST_MESH.DDy_LST=DDy_LST;
end

