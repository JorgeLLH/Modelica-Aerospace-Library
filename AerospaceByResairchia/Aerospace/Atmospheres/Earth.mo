within Aerospace.Atmospheres;

package Earth
  extends Aerospace.Icons.Packages;

  model EarthAtmosphere
    extends Icons.Models;
    parameter Modelica.Units.SI.Acceleration g0 = 9.80665;//Sea-level gravity
    parameter Real R(unit="J/(mol.K)") = Modelica.Constants.R * 10 ^ 3;//Gas constant 8.31432 * 10 ^ 3 -> Pressure in [Pa] thus the *10^3
    parameter Real beta(unit="kg/(s.m.K^1/2)") = 0.000001458;
    parameter Modelica.Units.SI.Temperature S = 110.4;//Sutherland constant
    type GeopotentialHeight = Real(unit = "m'");
    parameter GeopotentialHeight Hb[:] = {0, 11000, 20000, 32000, 47000, 51000, 71000, 84852};//Geopotential height atmospheric layers
    type MolecularScaleTemperatureGradient = Real(unit = "K/m'");
    parameter MolecularScaleTemperatureGradient L_Mb[:] = {-0.0065, 0.0, 0.001, 0.0028, 0.0, -0.0028, -0.0020};//Molecular scale temperature gradient
    parameter Real gamma = 1.4;//Specific heat ratio
    parameter Real M0(unit="kg/(kmol)") = 28.9644;//Molecular weight
    parameter Modelica.Units.SI.Length r0 = 6356766;//Earth radius
    parameter Modelica.Units.SI.Temperature Tb0 = 288.15;//Sea-level temperature
    parameter Modelica.Units.SI.Density rho0 = 1.224813066;//Sea-level density
    parameter Modelica.Units.SI.Velocity a0 = 340.3199579;//Sea-level speed of sound
    parameter Modelica.Units.SI.DynamicViscosity mu0 = 0.000174;//Sea-level dynamic viscosity
    Modelica.Units.SI.Length h;//Geopotential height
    Modelica.Units.SI.Length z;//Geometric height
    Modelica.Units.SI.Temperature Tb[7];//Layer temperature
    Modelica.Units.SI.Pressure Pb[7];//Layer pressure
    Integer b;//Atmospheric layer
    Modelica.Units.SI.Pressure P;//Local pressure
    Modelica.Units.SI.Temperature T;//Local temperature
    Modelica.Units.SI.Density rho;//Local density
    Modelica.Units.SI.Velocity a;//Local speed of sound
    Modelica.Units.SI.DynamicViscosity mu;//Local dynamic viscosity
    Modelica.Units.SI.SpecificHeatCapacity air_cp;//Local SpecificHeatCapacity
    Modelica.Units.SI.SpecificHeatCapacity air_cv;//Local SpecificHeatCapacity
    Modelica.Units.SI.SpecificEnthalpy air_h;//Local SpecificEnthalpy
    Modelica.Units.SI.SpecificEntropy air_s;//Local SpecificEntropy
    Modelica.Units.SI.RelativePressureCoefficient air_beta;//Local RelativePressureCoefficient
    Modelica.Units.SI.IsothermalCompressibility air_kappa;//Local IsothermalCompressibility
  equation
    z = 0;
    //Geopotential altitude correction
    h = r0 * z / (r0 + z);
    //Reference values computation
    Tb[1] = 288.15;
    for n in 1:6 loop
      Tb[n + 1] = Tb[n] + L_Mb[n] * (Hb[n + 1] - Hb[n]);
    end for;
    Pb[1] = 101325;
    for n in 1:6 loop
      if L_Mb[n] == 0 then
        Pb[n + 1] = Pb[n] * exp((-g0) * M0 * (Hb[n + 1] - Hb[n]) / (R * Tb[n]));
      else
        Pb[n + 1] = Pb[n] * (Tb[n] / (Tb[n] + L_Mb[n] * (Hb[n + 1] - Hb[n]))) ^ (g0 * M0 / (R * L_Mb[n]));
      end if;
    end for;
  //P and T values computation
    if h < Hb[2] then
      b = 1;
    elseif h < Hb[3] then
      b = 2;
    elseif h < Hb[4] then
      b = 3;
    elseif h < Hb[5] then
      b = 4;
    elseif h < Hb[6] then
      b = 5;
    elseif h < Hb[7] then
      b = 6;
    else
      b = 7;
    end if;
    if L_Mb[b] == 0 then
      P = Pb[b] * exp((-g0 * M0 * (h - Hb[b])) / (R * Tb[b]));
    else
      P = Pb[b] * (Tb[b] / (Tb[b] + L_Mb[b] * (h - Hb[b]))) ^ (g0 * M0 / (R * L_Mb[b]));
    end if;
    T = Tb[b] + L_Mb[b] * (h - Hb[b]);
    rho = P * M0 / (R * T);
    a = sqrt(gamma * R * T / M0);
    mu = beta * T ^ (3. / 2) / (T + S);
  
  //-------- Additional properties --------------------------
  
    air_cp = Modelica.Media.Air.ReferenceAir.Air_Utilities.cp_pT(P, T);
    air_cv = Modelica.Media.Air.ReferenceAir.Air_Utilities.cv_pT(P, T);
    air_h = Modelica.Media.Air.ReferenceAir.Air_Utilities.h_pT(P, T);
    air_s = Modelica.Media.Air.ReferenceAir.Air_Utilities.s_pT(P, T);
    air_beta = Modelica.Media.Air.ReferenceAir.Air_Utilities.beta_pT(P, T);
    air_kappa = Modelica.Media.Air.ReferenceAir.Air_Utilities.kappa_pT(P, T);
  
  annotation(Diagram(coordinateSystem(extent = {{-20, 20}, {20, -20}})));

  end EarthAtmosphere;
end Earth;
