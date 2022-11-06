within Aerospace;
package AircraftDesign
    extends Aerospace.Icons.Packages;

    package FixedWing
      extends Aerospace.Icons.Packages;

      package UnderDevelopment "Models under development. Use not recommended yet."
        extends Icons.InternalUse;

        model LiftingLineTheory
          extends Icons.Models;
          //References used:
          //1) "Fundamentals of Aerodynamics" by John D. Anderson.
          //2) "Aircraft Design a Systems Engineering Approach" by Mohammad H. Sadraey.
          //-----------------------
          parameter Integer N = 10;
          parameter Real AR = 8;
          parameter Real lambda = 0.6;
          parameter Modelica.Units.SI.Area S = 25;
          //-----------------------
          parameter Modelica.Units.SI.Position b = sqrt(AR * S);
          parameter Modelica.Units.SI.Position rootChord = 1.5 * (1 + lambda) * MAC / (1 + lambda + lambda ^ 2);
          parameter Modelica.Units.SI.Position tipChord = lambda*rootChord;
          //-----------------------
          parameter Modelica.Units.SI.Position MAC = S / b;
          parameter Modelica.Units.NonSI.Angle_deg alpha_twist = -1;
          parameter Modelica.Units.NonSI.Angle_deg rootPitch = 2;
          parameter Modelica.Units.NonSI.Angle_deg alpha_zero = -1.5;
          parameter Real C_lalpha = 6.3;
          //-----------------------
          parameter Real theta[:] = linspace(Modelica.Constants.pi / (2 * N), Modelica.Constants.pi / 2, N);
          parameter Modelica.Units.NonSI.Angle_deg alpha[:] = linspace(rootPitch + alpha_twist, rootPitch, N);
          parameter Modelica.Units.SI.Position z[:] = b / 2 * cos(theta);
          parameter Modelica.Units.SI.Position meanC[:] = array(rootChord * (1 - (1 - lambda) * cos(theta[i])) for i in 1:N);
          parameter Real z_adim[:] = b / 2 * cos(theta) / (b / 2);
          //-----------------------
          Real C_Lwing;
          Real C_Llocal[N];
          Real A[N];
          Real delta;
          Real e;
          Real C_Dinduced;
        equation
//------------------
          for i in 1:N loop
            meanC[i] * C_lalpha / (4 * b) * (alpha[i] - alpha_zero) / (180 / Modelica.Constants.pi) * sin(theta[i]) = sum(A[j] * sin((2 * j - 1) * theta[i]) * (meanC[i] * C_lalpha / (4 * b) * (2 * j - 1) + sin(theta[i])) for j in 1:N);
            C_Llocal[i] * meanC[i] = 4 * b * sum(A[j] * sin((2 * j - 1) * theta[i]) for j in 1:N);
          end for;
//------------------
          C_Lwing = A[1] * Modelica.Constants.pi * AR;
          delta = sum(i * (A[i] / A[1]) ^ 2 for i in 2:N);
          e = 1 / (1 + delta);
          C_Dinduced = C_Lwing * C_Lwing / (Modelica.Constants.pi * e * AR);
        end LiftingLineTheory;

      model WingForces
      extends LiftingLineTheory;
      Real Lift_local[N];
      equation
  for i in 1:N loop
          C_Llocal[i]*0.5 * 1.225 * V_inf^2 = Lift_local[i];
        end for;
      end WingForces;
        annotation(
          Icon);
      end UnderDevelopment;
      annotation(
        Icon(graphics = {Polygon(fillColor = {180, 180, 180}, fillPattern = FillPattern.Solid, points = {{-80, 20}, {-80, -20}, {0, -40}, {80, -20}, {80, 20}, {0, 40}, {-80, 20}}), Line(origin = {-60, 60}, points = {{0, 20}, {0, -20}, {0, -20}}), Line(origin = {60, 60}, points = {{0, 20}, {0, -20}}), Line(origin = {-19.9104, 59.9701}, points = {{0, 20}, {0, -20}, {0, -20}}), Line(origin = {19.8657, 59.9402}, points = {{0, 20}, {0, -20}, {0, -20}}), Line(origin = {-65, 45}, points = {{5, -5}, {5, -5}, {-5, 5}}), Line(origin = {-24.9104, 45.2836}, points = {{5, -5}, {5, -5}, {-5, 5}}), Line(origin = {14.8657, 45.2537}, points = {{5, -5}, {5, -5}, {-5, 5}}), Line(origin = {54.9553, 45.2238}, points = {{5, -5}, {5, -5}, {-5, 5}}), Line(origin = {-54.6866, 45.2836}, rotation = -90, points = {{5, -5}, {5, -5}, {-5, 5}}), Line(origin = {-14.597, 45.5672}, rotation = -90, points = {{5, -5}, {5, -5}, {-5, 5}}), Line(origin = {24.8657, 45.5373}, rotation = -90, points = {{5, -5}, {5, -5}, {-5, 5}}), Line(origin = {65.2687, 45.5074}, rotation = -90, points = {{5, -5}, {5, -5}, {-5, 5}})}));
    end FixedWing;

    package RotaryWing   "Models under development. Use not recommended yet."
      extends Aerospace.Icons.Packages;

      package UnderDevelopment
        extends Icons.InternalUse;

        model BEMT
          extends Icons.Models;
          //Parameters
          parameter Integer E(min = 1) = 50;
          parameter Real C_lalpha = 5.95;
          //2*2*Modelica.Constants.pi;
          parameter Real R = 0.9144;
          parameter Real R_cut_out = 0.1524 / 0.9144;
          parameter Real chord = 0.048006;
          parameter Real Blade_Pitch = 10.34 * Modelica.Constants.pi / 180;
          parameter Real Blade_Linear_Twist = -0 * Modelica.Constants.pi / 180;
          parameter Real r[:] = linspace(R_cut_out + delta_r / 2, 1 - delta_r / 2, E);
          parameter Real N_b = 3;
          parameter Real lambda_c = 0;
          parameter Real solidity = N_b * chord / (R * Modelica.Constants.pi);
          parameter Real delta_r = (1 - R_cut_out) / E;
          parameter Real PrandtlTipLoss_max_iter = 10;
          Real F[E];
          Real f[E];
          //Variables
          Real lambda_ideal[E];
          Real lambda[E];
          Real Ct_ideal;
          Real Ct;
          Integer iter_num;
        algorithm
          F := ones(E);
          f := ones(E);
          for i in 1:E loop
            lambda_ideal[i] := ((solidity * C_lalpha / (16 * F[i]) - lambda_c / 2) ^ 2 + solidity * C_lalpha * (Blade_Pitch + Blade_Linear_Twist * (r[i] - R_cut_out) / (1 - R_cut_out)) * r[i] / (8 * F[i])) ^ 0.5 - (solidity * C_lalpha / (16 * F[i]) - lambda_c / 2);
          end for;
          Ct_ideal := sum(solidity * (C_lalpha / 2) * ((Blade_Pitch + Blade_Linear_Twist * (r[i] - R_cut_out) / (1 - R_cut_out)) * r[i] * r[i] - lambda_ideal[i] * r[i]) * delta_r for i in 1:E);
//
          iter_num := 0;
          while F[E] > 0.001 and iter_num < PrandtlTipLoss_max_iter loop
            iter_num := iter_num + 1;
            for i in 1:E loop
              lambda[i] := ((solidity * C_lalpha / (16 * F[i]) - lambda_c / 2) ^ 2 + solidity * C_lalpha * (Blade_Pitch + Blade_Linear_Twist * (r[i] - R_cut_out) / (1 - R_cut_out)) * r[i] / (8 * F[i])) ^ 0.5 - (solidity * C_lalpha / (16 * F[i]) - lambda_c / 2);
              f[i] := N_b / 2 * ((1 - r[i]) / lambda[i]);
              F[i] := 2 / Modelica.Constants.pi * acos(exp(-f[i]));
            end for;
          end while;
        equation
          Ct = sum(solidity * (C_lalpha / 2) * ((Blade_Pitch + Blade_Linear_Twist * (r[i] - R_cut_out) / (1 - R_cut_out)) * r[i] * r[i] - lambda[i] * r[i]) * delta_r for i in 1:E);
        end BEMT;
      end UnderDevelopment;
      annotation(
        Icon(graphics = {Ellipse(extent = {{-80, 80}, {80, -80}}, endAngle = 360), Rectangle(origin = {50, 0}, fillColor = {180, 180, 180}, fillPattern = FillPattern.Solid, extent = {{-30, 8}, {30, -8}}), Rectangle(origin = {-50, 0}, fillColor = {180, 180, 180}, fillPattern = FillPattern.Solid, extent = {{-30, 8}, {30, -8}}), Rectangle(origin = {0, 50}, rotation = -90, fillColor = {180, 180, 180}, fillPattern = FillPattern.Solid, extent = {{-30, 8}, {30, -8}}), Rectangle(origin = {0, -50}, rotation = -90, fillColor = {180, 180, 180}, fillPattern = FillPattern.Solid, extent = {{-30, 8}, {30, -8}}), Line(points = {{-20, 0}, {20, 0}}), Line(points = {{0, 20}, {0, -20}})}));
    end RotaryWing;
    annotation(
      Icon(graphics = {Ellipse(origin = {0, 10}, fillColor = {180, 180, 180}, fillPattern = FillPattern.Solid, extent = {{-16, 70}, {16, -70}}), Polygon(origin = {0, -50}, fillColor = {180, 180, 180}, fillPattern = FillPattern.Solid, points = {{-12, 12}, {-30, -12}, {30, -12}, {12, 12}, {-12, 12}}), Rectangle(origin = {1, 25}, fillColor = {180, 180, 180}, fillPattern = FillPattern.Solid, extent = {{-67, 15}, {67, -15}})}));
end AircraftDesign;
