package Aerospace
// --------------------------------------------------------------------
  //                       Aerospace for Modelica
  //         Library for Aerospace Design for Modelica
// --------------------------------------------------------------------
// Copyright (C) 2022  Jorge L. Lovaco
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this software.  If not, see <http://www.gnu.org/licenses/>.
//
// Primary author: Jorge L. Lovaco <media@resairchia.com>
// Developed by: Resairchia
// Version: 0.0.001
// This is a beta version and still under development. Use at your own risk.
  extends Aerospace.Icons.Packages;

  package AircraftDesign
    extends Aerospace.Icons.Packages;

    package FixedWing
      extends Aerospace.Icons.Packages;

      package UnderDevelopment
        extends Icons.InternalUse;

        model LiftingLineTheory
        //References used:
          //1) "Fundamentals of Aerodynamics" by John D. Anderson.
          //2) "Aircraft Design a Systems Engineering Approach" by Mohammad H. Sadraey.
          parameter Modelica.Units.SI.Position b = sqrt(AR * S);
          parameter Modelica.Units.SI.Position rootChord = 1.5 * (1 + lambda) * MAC / (1 + lambda + lambda ^ 2);
          parameter Modelica.Units.SI.Position tipChord = 1.3;
          parameter Modelica.Units.SI.Area S = 25;
          //-----------------------
          parameter Integer N = 10;
          parameter Real AR = 8;
          parameter Real lambda = 0.6;
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
        annotation(
          Icon);
      end UnderDevelopment;
      annotation(
        Icon(graphics = {Polygon(fillColor = {180, 180, 180}, fillPattern = FillPattern.Solid, points = {{-80, 20}, {-80, -20}, {0, -40}, {80, -20}, {80, 20}, {0, 40}, {-80, 20}}), Line(origin = {-60, 60}, points = {{0, 20}, {0, -20}, {0, -20}}), Line(origin = {60, 60}, points = {{0, 20}, {0, -20}}), Line(origin = {-19.9104, 59.9701}, points = {{0, 20}, {0, -20}, {0, -20}}), Line(origin = {19.8657, 59.9402}, points = {{0, 20}, {0, -20}, {0, -20}}), Line(origin = {-65, 45}, points = {{5, -5}, {5, -5}, {-5, 5}}), Line(origin = {-24.9104, 45.2836}, points = {{5, -5}, {5, -5}, {-5, 5}}), Line(origin = {14.8657, 45.2537}, points = {{5, -5}, {5, -5}, {-5, 5}}), Line(origin = {54.9553, 45.2238}, points = {{5, -5}, {5, -5}, {-5, 5}}), Line(origin = {-54.6866, 45.2836}, rotation = -90, points = {{5, -5}, {5, -5}, {-5, 5}}), Line(origin = {-14.597, 45.5672}, rotation = -90, points = {{5, -5}, {5, -5}, {-5, 5}}), Line(origin = {24.8657, 45.5373}, rotation = -90, points = {{5, -5}, {5, -5}, {-5, 5}}), Line(origin = {65.2687, 45.5074}, rotation = -90, points = {{5, -5}, {5, -5}, {-5, 5}})}));
    end FixedWing;

    package RotaryWing
      extends Aerospace.Icons.Packages;
      annotation(
        Icon(graphics = {Ellipse(extent = {{-80, 80}, {80, -80}}, endAngle = 360), Rectangle(origin = {50, 0}, fillColor = {180, 180, 180}, fillPattern = FillPattern.Solid, extent = {{-30, 8}, {30, -8}}), Rectangle(origin = {-50, 0}, fillColor = {180, 180, 180}, fillPattern = FillPattern.Solid, extent = {{-30, 8}, {30, -8}}), Rectangle(origin = {0, 50}, rotation = -90, fillColor = {180, 180, 180}, fillPattern = FillPattern.Solid, extent = {{-30, 8}, {30, -8}}), Rectangle(origin = {0, -50}, rotation = -90, fillColor = {180, 180, 180}, fillPattern = FillPattern.Solid, extent = {{-30, 8}, {30, -8}}), Line(points = {{-20, 0}, {20, 0}}), Line(points = {{0, 20}, {0, -20}})}));
    end RotaryWing;
    annotation(
      Icon(graphics = {Ellipse(origin = {0, 10}, fillColor = {180, 180, 180}, fillPattern = FillPattern.Solid, extent = {{-16, 70}, {16, -70}}), Polygon(origin = {0, -50}, fillColor = {180, 180, 180}, fillPattern = FillPattern.Solid, points = {{-12, 12}, {-30, -12}, {30, -12}, {12, 12}, {-12, 12}}), Rectangle(origin = {1, 25}, fillColor = {180, 180, 180}, fillPattern = FillPattern.Solid, extent = {{-67, 15}, {67, -15}})}));
  end AircraftDesign;

  package Icons
  extends Aerospace.Icons.Packages;
    partial package Background
      annotation(
        Icon(graphics = {Rectangle(lineColor = {171, 171, 171},fillColor = {236, 236, 236}, fillPattern = FillPattern.VerticalCylinder,lineThickness = 1, extent = {{-100, 100}, {100, -100}})}));
    end Background;
  
    partial package Packages
    extends Background;
      annotation(
        Icon(graphics = {Polygon(origin = {65, -57}, lineColor = {85, 0, 255}, fillColor = {0, 255, 0}, fillPattern = FillPattern.Solid, lineThickness = 0.75, points = {{-17, -27}, {-17, 27}, {17, 27}, {17, 15}, {5, 15}, {5, 21}, {-7, 21}, {-7, -27}, {-17, -27}}), Polygon(origin = {83, -69}, lineColor = {85, 0, 255}, fillColor = {255, 170, 0}, fillPattern = FillPattern.Solid, lineThickness = 0.5, points = {{7, 21}, {1, 11}, {-1, -3}, {1, -15}, {-7, -27}, {-11, -27}, {-1, -17}, {-3, 1}, {1, 17}, {11, 25}, {13, 23}, {7, 21}}, smooth = Smooth.Bezier), Polygon(origin = {75, -69}, lineColor = {85, 0, 255}, fillColor = {255, 170, 0}, fillPattern = FillPattern.Solid, lineThickness = 0.5, points = {{1, 5}, {-3, 3}, {-3, -7}, {7, -11}, {5, -7}, {-1, -5}, {-1, 3}, {7, 5}, {7, 7}, {1, 5}}, smooth = Smooth.Bezier), Polygon(origin = {87, -73}, rotation = 180, lineColor = {85, 0, 255}, fillColor = {255, 170, 0}, fillPattern = FillPattern.Solid, lineThickness = 0.5, points = {{1, 5}, {-3, 3}, {-3, -7}, {7, -11}, {5, -7}, {-1, -5}, {-1, 3}, {7, 5}, {7, 7}, {1, 5}}, smooth = Smooth.Bezier)}));
    end Packages;
  
    partial package ExamplesPackageIcon
    extends Aerospace.Icons.Packages;
      annotation(
        Icon(graphics = {Polygon(origin = {-10, 0}, lineColor = {85, 0, 255}, fillColor = {0, 255, 0}, fillPattern = FillPattern.Solid, lineThickness = 2, points = {{-70, 60}, {70, 0}, {-70, -60}, {-70, 60}})}));
    end ExamplesPackageIcon;
  
    partial model ExampleModelIcon
      annotation(
        Icon(graphics = {Ellipse(lineColor = {255, 170, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 2, extent = {{-80, 80}, {80, -80}}), Polygon(origin = {12, 0}, lineColor = {85, 0, 255}, fillColor = {0, 255, 0}, fillPattern = FillPattern.Solid, lineThickness = 2, points = {{-70, 60}, {70, 0}, {-70, -60}, {-70, 60}})}));
    end ExampleModelIcon;
  
    partial package InternalUse
      extends Background;
    annotation(
        Icon(graphics = {Ellipse(fillPattern = FillPattern.Solid,lineThickness = 2, extent = {{-80, 80}, {80, -80}}), Ellipse(fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid,lineThickness = 2, extent = {{-60, 60}, {60, -60}}), Polygon(origin = {-2, -2}, fillPattern = FillPattern.Solid, points = {{-60, 40}, {-40, 60}, {60, -40}, {40, -60}, {-60, 40}, {-60, 40}})}));
    end InternalUse;
    annotation(
      Icon(coordinateSystem(initialScale = 0.1), graphics = {Text(extent = {{-80, -80}, {80, 80}}, textString = "Ω")}));
  end Icons;

  package LandingGearDesign
    extends Aerospace.Icons.Packages;

    package Tyres
  extends Aerospace.Icons.Packages;
    end Tyres;

    package Skids
  extends Icons.Packages;
    end Skids;
  end LandingGearDesign;
  annotation(
    Icon(graphics = {Text(origin = {0, 10}, extent = {{-66, 84}, {66, -84}}, textString = "A")}),
    uses(Modelica(version = "4.0.0")));
end Aerospace;
