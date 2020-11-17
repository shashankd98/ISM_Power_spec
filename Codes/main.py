import TB_vs_vlsr
import Col_density_plot
import Multiphase_ISM
import Planck
import Solenoidal_fraction
import WISE

print("\nPlanck Survey\n")
Planck.main()
print("\nWISE Survey\n")
WISE.main()
print("\nLAB Survey\n")
TB_vs_vlsr.main()
print("\nLAB Survey-Column Density Plots\n")
Col_density_plot.main()
print("\nLAB Survey-Multiphase ISM\n")
Multiphase_ISM.main()
print("\nOrion B Survey-Solenoidal Fraction")
Solenoidal_fraction.main()