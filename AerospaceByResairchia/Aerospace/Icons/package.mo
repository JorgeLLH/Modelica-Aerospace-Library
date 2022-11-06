// --------------------------------------------------------------------
// Resairchia Icons for Libraries in Modelica
// --------------------------------------------------------------------
// Copyright (C) 2021  Jorge L. Lovaco
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
within Aerospace;

package Icons
    //                       Icons for Modelica
                                      extends Packages;

  partial package Background
    annotation(
      Icon(graphics = {Rectangle(lineColor = {171, 171, 171}, fillColor = {236, 236, 236}, fillPattern = FillPattern.VerticalCylinder, lineThickness = 1, extent = {{-100, 100}, {100, -100}})}));
  end Background;

  partial package Packages
    extends Background;
    annotation(
      Icon(graphics = {Polygon(origin = {65, -57}, lineColor = {85, 0, 255}, fillColor = {0, 255, 0}, fillPattern = FillPattern.Solid, lineThickness = 0.75, points = {{-17, -27}, {-17, 27}, {17, 27}, {17, 15}, {5, 15}, {5, 21}, {-7, 21}, {-7, -27}, {-17, -27}}), Polygon(origin = {83, -69}, lineColor = {85, 0, 255}, fillColor = {255, 170, 0}, fillPattern = FillPattern.Solid, lineThickness = 0.5, points = {{7, 21}, {1, 11}, {-1, -3}, {1, -15}, {-7, -27}, {-11, -27}, {-1, -17}, {-3, 1}, {1, 17}, {11, 25}, {13, 23}, {7, 21}}, smooth = Smooth.Bezier), Polygon(origin = {75, -69}, lineColor = {85, 0, 255}, fillColor = {255, 170, 0}, fillPattern = FillPattern.Solid, lineThickness = 0.5, points = {{1, 5}, {-3, 3}, {-3, -7}, {7, -11}, {5, -7}, {-1, -5}, {-1, 3}, {7, 5}, {7, 7}, {1, 5}}, smooth = Smooth.Bezier), Polygon(origin = {87, -73}, rotation = 180, lineColor = {85, 0, 255}, fillColor = {255, 170, 0}, fillPattern = FillPattern.Solid, lineThickness = 0.5, points = {{1, 5}, {-3, 3}, {-3, -7}, {7, -11}, {5, -7}, {-1, -5}, {-1, 3}, {7, 5}, {7, 7}, {1, 5}}, smooth = Smooth.Bezier)}));
  end Packages;

  partial package ExamplesPackageIcon
    extends Icons.Packages;
    annotation(
      Icon(graphics = {Polygon(origin = {-10, 0}, lineColor = {85, 0, 255}, fillColor = {0, 255, 0}, fillPattern = FillPattern.Solid, lineThickness = 2, points = {{-70, 60}, {70, 0}, {-70, -60}, {-70, 60}})}));
  end ExamplesPackageIcon;

  partial model ExampleModelIcon
    annotation(
      Icon(graphics = {Ellipse(lineColor = {255, 170, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 2, extent = {{-80, 80}, {80, -80}}), Polygon(origin = {12, 0}, lineColor = {85, 0, 255}, fillColor = {0, 255, 0}, fillPattern = FillPattern.Solid, lineThickness = 2, points = {{-70, 60}, {70, 0}, {-70, -60}, {-70, 60}})}));
  end ExampleModelIcon;

  partial package InternalUse
    extends Icons.Background;
    annotation(
      Icon(graphics = {Ellipse(fillPattern = FillPattern.Solid, lineThickness = 2, extent = {{-80, 80}, {80, -80}}), Ellipse(fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 2, extent = {{-60, 60}, {60, -60}}), Polygon(origin = {-2, -2}, fillPattern = FillPattern.Solid, points = {{-60, 40}, {-40, 60}, {60, -40}, {40, -60}, {-60, 40}, {-60, 40}})}));
  end InternalUse;

  partial function functionIcon
    annotation(
      Icon(graphics = {Ellipse(lineColor = {85, 0, 255}, fillColor = {255, 170, 0}, fillPattern = FillPattern.Solid, lineThickness = 2, extent = {{-60, 60}, {60, -60}}), Text(origin = {-3, 5}, extent = {{-33, 31}, {33, -31}}, textString = "f")}));
  end functionIcon;

partial function derivativeIcon
  annotation(
    Icon(coordinateSystem(initialScale = 0.1), graphics = {Text(origin = {-41, 54}, extent = {{-33, 24}, {101, -124}}, textString = "D"), Text(origin = {24, -14}, extent = {{-22, 14}, {54, -38}}, textString = "x")}));
end derivativeIcon;

partial function integralIcon
  annotation(
    Icon(coordinateSystem(initialScale = 0.1), graphics = {Text(origin = {-43, 70}, extent = {{-31, 18}, {117, -130}}, textString = "∫")}),
    Diagram);
end integralIcon;

  annotation(
    Icon(coordinateSystem(initialScale = 0.1), graphics = {Text(extent = {{-80, -80}, {80, 80}}, textString = "Ω")}));
end Icons;
