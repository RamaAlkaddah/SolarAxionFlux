



# changes to the code: 

### In constants.hpp:

```
		const double m12_keV = 12.0e3; // 12 keV in eV
		const double m120_keV = 120.0e3; // 120 keV in eV
		const double m0_keV = 0.0e3; // 0 keV in eV
```

### In Solar_model.hpp
```
// Struct for resonant transverse plasmon production rates at different axion masses

struct GammaTPResonant {

  double m120;    // Gamma_TP for m_a = 120 eV

  double m12;     // Gamma_TP for m_a = 12 eV

  double m0;      // Gamma_TP for m_a = 0 eV

};
```

also around line 96: 

```
    // Resonant (massive) transverse plasmon production rates

    GammaTPResonant Gamma_TP_resonant(double omega, double r);

    double Gamma_TP_resonant_m120(double omega, double r);

    double Gamma_TP_resonant_m12(double omega, double r);

    double Gamma_TP_resonant_m0(double omega, double r);
```

### In test.hpp
inside void , line 21

```
 // Set true to run only TP + resonant TP calculations.

  const bool run_tp_only = true;
  
```

after t4s:
```
  

  // NEW: Resonant (massive) TP calculations using Gamma_TP_resonant (AGSS09 TP)

  auto t_m120_s = time_now();

  std::cout << "\n# Calculating resonant transversal plasmon spectrum (m_a = 120 eV)..." << std::endl;

  fully_integrate_d2Phi_a_domega_drho_in_rho(test_ergs, s, &SolarModel::Gamma_TP_resonant_m120, output_path + "TP_resonant_m120.dat");

  auto t_m120_e = time_now();

  std::cout << "# Calculating resonant TP (120 eV) took " << duration_cast<seconds>(t_m120_e-t_m120_s).count() << " seconds." << std::endl;

  

  auto t_m12_s = time_now();

  std::cout << "\n# Calculating resonant transversal plasmon spectrum (m_a = 12 eV)..." << std::endl;

  fully_integrate_d2Phi_a_domega_drho_in_rho(test_ergs, s, &SolarModel::Gamma_TP_resonant_m12, output_path + "TP_resonant_m12.dat");

  auto t_m12_e = time_now();

  std::cout << "# Calculating resonant TP (12 eV) took " << duration_cast<seconds>(t_m12_e-t_m12_s).count() << " seconds." << std::endl;

  

  auto t_m0_s = time_now();

  std::cout << "\n# Calculating resonant transversal plasmon spectrum (m_a = 0 eV)..." << std::endl;

  fully_integrate_d2Phi_a_domega_drho_in_rho(test_ergs, s, &SolarModel::Gamma_TP_resonant_m0, output_path + "TP_resonant_m0.dat");

  auto t_m0_e = time_now();

  std::cout << "# Calculating resonant TP (0 eV) took " << duration_cast<seconds>(t_m0_e-t_m0_s).count() << " seconds." << std::endl;
```

## Inside results / src 
### In python_wrapper.cpp:
```
 .def("tp_resonant_m120", pybind11::vectorize(&SolarModel::Gamma_TP_resonant_m120), "Resonant transverse plasmon production rate (m_a = 120 eV)", "omega"_a, "radius"_a)

    .def("tp_resonant_m12", pybind11::vectorize(&SolarModel::Gamma_TP_resonant_m12), "Resonant transverse plasmon production rate (m_a = 12 eV)", "omega"_a, "radius"_a)

    .def("tp_resonant_m0", pybind11::vectorize(&SolarModel::Gamma_TP_resonant_m0), "Resonant transverse plasmon production rate (m_a = 0 eV)", "omega"_a, "radius"_a)
```

### In solar_model.cpp 
after gamma_TP
```
GammaTPResonant SolarModel::Gamma_TP_resonant(double omega, double r) {

  static const double geom_factor = 1.0;

  static const double photon_polarization = 2.0;

  

  // --- fixed axion masses (in eV) ---

  const double m_eV_120 = 0.120;

  const double m_eV_12  = 0.012;

  const double m_eV_0   = 0.000;

  

  GammaTPResonant out{0.00, 0.00, 0.00};

  

  double om2   = omega * omega;

  double om_pl_sq = omega_pl_squared(r);

  


  // no TP mode below plasma frequency

  if (om_pl_sq > om2) {

    return out;

  }

  

  double T = temperature_in_keV(r);

  double u = omega / T;

  

  // photon damping rate from opacity (frequency-dependent)

  double gamma = -gsl_expm1(-u) * opacity(omega, r);

  

  // average transverse B^2

  double average_b_field_sq = gsl_pow_2(bfield(r)) / 3.0;

  

  // Δ_T^2 = g_{aγ}^2 ⟨B_T^2⟩ / 4

  double DeltaTsq = g_agg * g_agg * average_b_field_sq / 4.0;

  

  // Δ_γ = -ω_pl^2 / (2ω)

  double Delta_gamma = -om_pl_sq / (2.0 * omega);

  

  auto rate_for_mass = [&](double m_keV) -> double {

    double m2 = m_keV * m_keV;

  

    // Δ_a = -m_a^2 / (2ω)

    double Delta_a = -m2 / (2.0 * omega);

  

    // Δ_γ - Δ_a

    double Delta_diff = Delta_gamma - Delta_a;

  

    // Breit–Wigner denominator

    double denom = gsl_pow_2(Delta_diff) + gsl_pow_2(0.5 * gamma);

  

    // full TP rate (same thermal factor as your code)

    double result = geom_factor * photon_polarization

                    * gamma * DeltaTsq

                    / (denom * gsl_expm1(u));

  

    return result;

  };

  

  out.m120 = rate_for_mass(m_eV_120);

  out.m12  = rate_for_mass(m_eV_12);

  out.m0   = rate_for_mass(m_eV_0);

  

  return out;

}

  

// Wrapper functions to extract individual masses from Gamma_TP_resonant for integration

double SolarModel::Gamma_TP_resonant_m120(double omega, double r) {

  return Gamma_TP_resonant(omega, r).m120;

}

  

double SolarModel::Gamma_TP_resonant_m12(double omega, double r) {

  return Gamma_TP_resonant(omega, r).m12;

}

  

double SolarModel::Gamma_TP_resonant_m0(double omega, double r) {

  return Gamma_TP_resonant(omega, r).m0;

}
```
## in plot_results.py 

in the beginning
```
    # Load resonant TP files (Breit-Wigner resonant, AGSS09 TP) - with fallback for missing files

    try:

        TP_resonant_m120 = np.genfromtxt(script_path+"/TP_resonant_m120.dat")

    except:

        TP_resonant_m120 = None

        print("WARNING: TP_resonant_m120.dat not found, skipping this file.")

    try:

        TP_resonant_m12 = np.genfromtxt(script_path+"/TP_resonant_m12.dat")

    except:

        TP_resonant_m12 = None

        print("WARNING: TP_resonant_m12.dat not found, skipping this file.")

    try:

        TP_resonant_m0 = np.genfromtxt(script_path+"/TP_resonant_m0.dat")

    except:

        TP_resonant_m0 = None

        print("WARNING: TP_resonant_m0.dat not found, skipping this file.")

except OSError as e:

    print("\nResult files from standard tests not found! You need to run the 'test_library'"

           "program from the 'bin/' folder first.\n")

    raise
```

then to plot 

```
  

## Massive TP comparison plot (AGSS09 TP)

fig, ax = plt.subplots()

plot_setup()

ax.plot(res6[:,0], res6[:,1]/1.0e10, 'k--', label=r'TP (m$_a$ = 0 eV, Off-Res)')

if TP_resonant_m120 is not None:

    ax.plot(TP_resonant_m120[:,0], TP_resonant_m120[:,1]/1.0e10, '-', color='#1f77b4', linewidth=2,alpha=0.7, label=r'TP (m$_a$ = 120 eV)')

if TP_resonant_m12 is not None:

    ax.plot(TP_resonant_m12[:,0], TP_resonant_m12[:,1]/1.0e10, '-', color="#eb327f", linewidth=2,alpha=0.7, label=r'TP (m$_a$ = 12 eV)')

if TP_resonant_m0 is not None:

    ax.plot(TP_resonant_m0[:,0], TP_resonant_m0[:,1]/1.0e10, '-', color='#2ca02c', linewidth=2,alpha=0.7, label=r'TP (m$_a$ = 0 eV)')

  

ax.set_title(r'Transverse plasmon mode: Effect of axion mass (AGSS09 TP Rosseland), $g_{a\gamma\gamma} = \SI{1e-10}{\GeV^{-1}}$')

ax.set_xlabel(r'Energy $\omega$ [keV]')

ax.set_ylabel(r'Axion flux $\mathrm{d}\Phi_a/\mathrm{d}\omega$ [\SI{e10}{\per\cm\squared\per\keV\per\s}]')

ax.set_xlim([0.1,15])

ax.set_ylim([0.01,1.0e4])

ax.set_yscale('log')

ax.set_xscale('log')

ax.legend(frameon=False, loc='best')

#ax.grid(True, alpha=0.3, which='both')

  

plt.savefig(script_path+"/massive_TP_comparison.pdf", bbox_inches='tight')

#plt.show()

plt.close()

plt.close()
```

