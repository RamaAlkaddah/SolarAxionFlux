
### 1. in constants.hpp: 
change the coupling constant to
	`const double g_agg = 0.5E-16`
### 2. in  solar_model.hpp: 
introduce these two function for the two masses in line 94
\\in public (the beginning of the code) add:
`struct GammaTPResonant {`

      `double m200;    // Gamma_TP for m_a = 200 eV`

      `double m12;     // Gamma_TP for m_a = 12 eV`

      `double m0;      // Gamma_TP for m_a = 0 eV`

    `};`


    double Gamma_TP_resonant_m131(double omega, double r);

    double Gamma_TP_resonant_m11(double omega, double r);

	Note: m0 (rsonant) turned out to be negligable so we can skip it 
### 3. in test.hpp: 
change the energy counts and range to: 
```
 

  const int n_erg_values = 700;

 ////

  std:: vector<double> test_ergs;

  for (int k=0; k<n_erg_values; k++) { test_ergs.push_back(0.001 + k*15.999/n_erg_values); }
```

and add test function for each mass:

```
  

  auto t_m131_s = time_now();

  std::cout << "\n# Calculating resonant transversal plasmon spectrum (m_a = 131 eV)..." << std::endl;

  fully_integrate_d2Phi_a_domega_drho_in_rho(test_ergs, s, &SolarModel::Gamma_TP_resonant_m131, output_path + "TP_resonant_m131.dat");

  auto t_m131_e = time_now();

  std::cout << "# Calculating resonant TP (131 eV) took " << duration_cast<seconds>(t_m131_e-t_m131_s).count() << " seconds." << std::endl;

  auto t_m11_s = time_now();

  std::cout << "\n# Calculating resonant transversal plasmon spectrum (m_a = 11 eV)..." << std::endl;

  fully_integrate_d2Phi_a_domega_drho_in_rho(test_ergs, s, &SolarModel::Gamma_TP_resonant_m11, output_path + "TP_resonant_m11.dat");

  auto t_m11_e = time_now();

  std::cout << "# Calculating resonant TP (11 eV) took " << duration_cast<seconds>(t_m11_e-t_m11_s).count() << " seconds." << std::endl;
```

### 4. in solar_moel.cpp: 
add a function that calculate the conversion rate for the massive TP then apply it on each mass in different function to avoid adding the mass as a new parameter (here i did each mass separetly to avoid some issues, i might come back to edit them later):
```


double SolarModel::Gamma_TP_resonant_m131(double omega, double r) {

    static const double geom_factor = 1.0;

    static const double photon_polarization = 2.0;

    const double m_keV = 0.131;  // 131 eV in keV

    // Full Breit-Wigner calculation inline

    if (omega < 1e-6) return 0;

    double T = temperature_in_keV(r);

    double u = omega / T;

    double gamma = -gsl_expm1(-u) * interpolate_rosseland_opacity(r);

     // --- NEW: narrow resonance window around omega_pl(r) = m_a ---

    gamma = std::max(gamma, 1e-4);

    double xi2 = gamma * omega;

    double om2 = omega * omega;

    double fwhm = (om2 > xi2) ? (sqrt(om2 + xi2) - sqrt(om2 - xi2)) : sqrt(om2 + xi2);

    double om_pl_sq = omega_pl_squared(r);

    if (std::abs(std::sqrt(om_pl_sq) - m_keV) > 18.0 * fwhm) { return 0; }

    // -----------------------------------------------------------

    double average_b_field_sq = gsl_pow_2(bfield(r)) / 3.0;

    double DeltaTsq = g_agg * g_agg * average_b_field_sq / 4.0;

    double Delta_gamma = -om_pl_sq / (2.0 * omega);

    double m2 = m_keV * m_keV;

    double Delta_a = -m2 / (2.0 * omega);

    double Delta_diff = Delta_gamma - Delta_a;

    double denom = gsl_pow_2(Delta_diff) + gsl_pow_2(0.5 * gamma);

    return geom_factor * photon_polarization * gamma * DeltaTsq / (denom * gsl_expm1(u));

}

  

double SolarModel::Gamma_TP_resonant_m11(double omega, double r) {

    static const double geom_factor = 1.0;

    static const double photon_polarization = 2.0;

    const double m_keV = 0.011;  // 11 eV in keV

    // Full Breit-Wigner calculation inline

    if (omega < 1e-6) return 0;

    double T = temperature_in_keV(r);

    double u = omega / T;

    double gamma = -gsl_expm1(-u) * interpolate_rosseland_opacity(r);

     // --- NEW: narrow resonance window around omega_pl(r) = m_a ---

    gamma = std::max(gamma, 1e-4);

    double xi2 = gamma * omega;

    double om2 = omega * omega;

    double fwhm = (om2 > xi2) ? (sqrt(om2 + xi2) - sqrt(om2 - xi2)) : sqrt(om2 + xi2);

    double om_pl_sq = omega_pl_squared(r);

    if (std::abs(std::sqrt(om_pl_sq) - m_keV) > 18.0 * fwhm) { return 0; }

    // -----------------------------------------------------------

    double average_b_field_sq = gsl_pow_2(bfield(r)) / 3.0;

    double DeltaTsq = g_agg * g_agg * average_b_field_sq / 4.0;

    double Delta_gamma = -om_pl_sq / (2.0 * omega);

    double m2 = m_keV * m_keV;

    double Delta_a = -m2 / (2.0 * omega);

    double Delta_diff = Delta_gamma - Delta_a;

    double denom = gsl_pow_2(Delta_diff) + gsl_pow_2(0.5 * gamma);

    return geom_factor * photon_polarization * gamma * DeltaTsq / (denom * gsl_expm1(u));

}
```

### 5. in pyhton_wrapper.cpp: 
add
```
    .def("tp_resonant_m131", pybind11::vectorize(&SolarModel::Gamma_TP_resonant_m131), "Resonant transverse plasmon production rate (m_a = 131 eV)", "omega"_a, "radius"_a)

    .def("tp_resonant_m11", pybind11::vectorize(&SolarModel::Gamma_TP_resonant_m11), "Resonant transverse plasmon production rate (m_a = 11 eV)", "omega"_a, "radius"_a)
```

---


### 6. For plotting: plot_results.py  
#### - for plotting TP 
```

## Massive TP comparison plot (AGSS09 TP)

# ==========================================

# Magnetic-field uncertainty helper, Flux scales as B^2

def magnetic_field_band(flux):

    B_low  = 0.25 * flux    # low  = (0.5B)^2 = 0.25

    B_high = 2.25 * flux   # high = (1.5B)^2 = 2.25

    return B_low, B_high

  

fig, ax = plt.subplots()

plot_setup()

# TP (m_a = 0 eV)

m0 = (res8[:,1] / scale)  #this is  off resonant one based on Rosseland, but we want to show the magnetic field uncertainty band

B_low, B_high = magnetic_field_band(m0)

ax.fill_between(

    res8[:,0],

    B_low,

    B_high,

    color="#4EAF21",

    alpha=0.20,

    linewidth=0

)

ax.plot(

    res8[:,0],

    m0,

    'k--',

    color="#000000",

    linewidth=2,

    alpha=0.7, label=r'TP (m$_a$ = 0 eV)'

)

# =========================

# TP (m_a = 131 eV)

m131 = (TP_resonant_m131[:,1] / scale)

B_low, B_high = magnetic_field_band(m131)

ax.fill_between(

    TP_resonant_m131[:,0],

    B_low,

    B_high,

    color='#1f77b4',

    alpha=0.20,

    linewidth=0

)

ax.plot(

    TP_resonant_m131[:,0],

    m131,

    '-',

    color='#1f77b4',

    linewidth=2,

    alpha=0.7, label=r'TP (m$_a$ = 131 eV)'

)

# =========================

# TP (m_a = 11 eV)

m11 = (TP_resonant_m11[:,1] / scale)

B_low, B_high = magnetic_field_band(m11)

ax.fill_between(

    TP_resonant_m11[:,0],

    B_low,

    B_high,

    color="#eb327f",

    alpha=0.20,

    linewidth=0

)

ax.plot(

    TP_resonant_m11[:,0],

    m11,

    '-',

    color="#eb327f",

    linewidth=2,

    alpha=0.7, label=r'TP (m$_a$ = 11 eV)'

)

  

ax.set_xlabel(r'$\omega$ [keV]')

ax.set_ylabel(r'$\mathrm{d}\Phi_a/\mathrm{d}\omega$ ' r'[\SI{e10}{\per\cm\squared\per\keV\per\s}]')

ax.set_xlim([0.001, 20])

ax.set_ylim([0.001, 2e3])

ax.set_yscale('log')

ax.set_xscale('log')

  

ax.legend(frameon=False, loc='upper left')

  

# ax.grid(True, alpha=0.3, which='both')

  

plt.savefig(script_path + "/massive_TP_comparison.pdf", bbox_inches='tight')

  

plt.show()

  

plt.close()

#plt.close()
```
