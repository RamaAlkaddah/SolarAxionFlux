# Comprehensive Comparison: Giannotti Non-Resonant TP vs. Hoof Rosseland Model

## 1. PHYSICS BACKGROUND

### 1.1 Physical Process: Coherent Conversion of Thermal Photons

Both Giannotti et al. (PhysRevD.102.123024) and Hoof et al. (arXiv:2101.08789) describe the same fundamental physics: the **coherent conversion of thermal photons to axionlike particles (ALPs) in large-scale solar magnetic fields**.

**Key Physical Features:**
- **Non-resonant zero-mass regime** (both methods):
  - Applicable when axion mass $m_a \approx 0$ (or $m_a \ll \omega_p$)
  - Dominated by transverse plasmon (TP) production
  - Coherent scattering in background magnetic field

- **Mechanism**: Thermal photon + magnetic field $\rightarrow$ axion
  - Photon virtual oscillation (energy-dependent damping rate)
  - Coherent momentum transfer via plasmon dispersion
  - Resonance broadened by thermal damping $\gamma$

### 1.2 Mathematical Framework

#### Giannotti Formulation
The Giannotti paper derives TP flux from coherent photon-axion conversion:

$$\Gamma_{\text{TP}}(\omega, r) = N_{\text{pol}} \cdot \frac{\Delta_T^2 \cdot \gamma(\omega,r)}{[\Delta_p(\omega) - \Delta_a(\omega)]^2 + [\gamma(\omega,r)/2]^2} \cdot \frac{1}{e^{\omega/T}-1}$$

Where:
- $N_{\text{pol}} = 2$ (two transverse polarization states)
- $\Delta_T^2 = \frac{g_{a\gamma}^2 \langle B_T^2 \rangle}{4}$ (transverse field coupling squared)
- $\gamma(\omega,r) = -\text{expm1}(-\omega/T) \cdot \kappa(\omega,r)$ (photon damping rate)
- $\Delta_p(\omega) = -\frac{\omega_p^2}{2\omega}$ (photon self-energy in plasma)
- $\Delta_a(\omega) = -\frac{m_a^2}{2\omega}$ (axion self-energy)
- For $m_a = 0$: $\Delta_a = 0$ (off-resonance condition)

#### Hoof Implementation
The Hoof code implements the non-resonant TP production via two approaches:

**A) Standard OP (Opacity Project) approach:**
```
Gamma_TP(ω, r) = geom_factor × photon_polarization × γ(ω,r) × ΔT² 
                / [(ΔP² + (γ/2)²) × expm1(ω/T)]
```

Where:
- $\gamma = -\text{expm1}(-u) \cdot \kappa(\omega, r)$ with $u = \omega/T$
- $\Delta_P^2 = \omega^2 [(\sqrt{1-\omega_p^2/\omega^2}-1)]^2$ (momentum transfer squared)
- $\Delta_T^2 = \frac{g_{a\gamma}^2 \langle B^2 \rangle/3}{4}$ (isotropic averaging)
- $\kappa(\omega,r)$ = frequency-dependent opacity from OP tables

**B) Rosseland Mean Opacity approach:**
```
Gamma_TP_Rosseland(ω, r) = geom_factor × photon_polarization × γ(ω,r) × ΔT²
                          / [(ΔP² + (γ/2)²) × expm1(ω/T)]
```

Identical to OP version, but:
- $\kappa(\omega,r) \rightarrow \kappa_R(r)$ (Rosseland mean opacity, frequency-independent)
- Single radius-dependent value for all photon energies

---

## 2. KEY PHYSICAL DIFFERENCES

### 2.1 Opacity Treatment

| Aspect | Giannotti | Hoof (OP) | Hoof (Rosseland) |
|--------|-----------|-----------|------------------|
| **Opacity Type** | Frequency-dependent | Frequency-dependent | Frequency-independent |
| **Energy Dependence** | $\kappa(\omega, r)$ | $\kappa(\omega, r)$ | $\kappa_R(r)$ only |
| **Temperature Dependence** | Implicit in $\kappa(\omega,T)$ | Implicit in OP tables | Explicit via Rosseland integral |
| **Physical Basis** | Electron-photon interactions | Thomson + bound-free absorption | Energy-weighted average opacity |
| **Valid Regime** | All photon energies | All photon energies | Radiation-dominated transport |

**Physical Implication:**
- **OP approach**: More accurate for individual frequency bands but computationally expensive; requires full energy-dependent interpolation
- **Rosseland approach**: Represents average opacity across spectrum; simplifies calculation but loses frequency resolution; particularly accurate for deep solar interior where radiation dominates

### 2.2 Momentum Transfer and Plasma Effects

| Feature | Giannotti | Hoof (OP/Rosseland) |
|---------|-----------|-------------------|
| **Momentum Transfer** | $\Delta_p = -\omega_p^2/(2\omega)$ | $\Delta_P = \omega[\sqrt{1-\omega_p^2/\omega^2}-1]$ |
| **Expression Type** | Simplified form (thin plasma limit) | Full relativistic form |
| **Plasma Wave Vector** | Implicit in dispersion | Explicit via $\omega_p$ |
| **Energy Threshold** | Soft ($\omega > \omega_p$ technically) | Hard ($\omega > \omega_p$ enforced) |

**Physics Comparison:**

The momentum transfer expressions differ significantly:
$$\text{Giannotti: } \Delta_p^2 = \frac{\omega_p^4}{4\omega^2}$$
$$\text{Hoof: } \Delta_P^2 = \omega^2 \left(\sqrt{1-\frac{\omega_p^2}{\omega^2}}-1\right)^2$$

For $\omega \gg \omega_p$ (high energy limit):
- Both $\sim \omega^2 \cdot O(1)$ 
- Agreement improves at high energies

For $\omega \sim \omega_p$ (resonance region):
- Hoof: $\Delta_P^2 \to 0$ rapidly
- Giannotti: More gradual approach
- **This difference affects low-energy flux significantly**

### 2.3 Magnetic Field Treatment

| Aspect | Giannotti | Hoof (Both OP/Rosseland) |
|--------|-----------|--------------------------|
| **Field Description** | Large-scale magnetic field | Parametrized dipole/toroidal field |
| **Averaging** | $\langle B_T^2 \rangle$ (unspecified) | $\langle B^2 \rangle = B^2/3$ (isotropic) |
| **Spatial Profile** | Coherent across domain | Stratified (radial, tachocline, outer) |
| **Values** | Calibrated to observations | Tunable parameters (default: $B_{rad}=200$ G, etc.) |

**Numerical Impact:**
- Hoof uses explicit magnetic field model with three zones
- Giannotti may use different field magnitude/profile
- This contributes to overall flux normalization differences

---

## 3. NUMERICAL IMPLEMENTATION COMPARISON

### 3.1 Key Differences in Code Structure

#### Hoof Implementation (from [src/solar_model.cpp](src/solar_model.cpp#L1115))

```cpp
double SolarModel::Gamma_TP(double omega, double r) {
  static const double geom_factor = 1.0;
  static const double photon_polarization = 2.0;
  
  if (omega_pl_squared(r) > omega*omega) { return 0; }
  
  double u = omega/temperature_in_keV(r);
  double gamma = -gsl_expm1(-u)*opacity(omega, r);
  double DeltaPsq = omega*omega * gsl_pow_2(
    sqrt(1-omega_pl_squared(r)/(omega*omega))-1
  );
  
  double average_b_field_sq = gsl_pow_2(bfield(r)) / 3.0;
  double DeltaTsq = g_agg*g_agg * average_b_field_sq / 4.0;
  
  double result = geom_factor * photon_polarization * gamma * DeltaTsq 
                / ( (DeltaPsq+gsl_pow_2(0.5*gamma)) * gsl_expm1(u) );
  return result;
}
```

**Key Features:**
1. Dimensionless units (energies in keV)
2. Early exit if $\omega < \omega_p$ (plasma frequency cutoff)
3. Explicit isotropic field averaging
4. Uses GSL special functions for numerical stability (`gsl_expm1`, `gsl_pow_2`)

#### Rosseland Variant (from [src/solar_model.cpp](src/solar_model.cpp#L1134))

```cpp
double SolarModel::Gamma_TP_Rosseland(double omega, double r) {
  // Identical structure to Gamma_TP()
  // Only difference: opacity source
  double gamma = -gsl_expm1(-u)*interpolate_rosseland_opacity(r);
  // ... rest identical
}
```

### 3.2 Numeric Stability Considerations

**Hoof's Stability Features:**
1. **`gsl_expm1(z)` instead of `exp(z)-1`**:
   - Prevents catastrophic cancellation for small $z = \omega/T$
   - Crucial for sub-keV photons where $z \ll 1$
   - Example: For $\omega = 0.1$ keV, $T = 1$ keV: $z = 0.1$
     - Direct: $\exp(0.1) - 1 = 1.10517 - 1 = 0.10517$ (loss of precision)
     - `expm1(0.1)` = 0.10517 (full precision)

2. **`gsl_pow_2` for squaring**:
   - Inline optimization for $x^2$ operations
   - Avoids overflow/underflow for intermediate values

3. **Early plasma frequency check**:
   - Prevents NaN from $\sqrt{\text{negative}}$ when $\omega < \omega_p$
   - Returns exactly 0 (no numerical artifacts)

### 3.3 Integration Strategy

**Spectral Integration (from [tests.hpp](include/solaxflux/tests.hpp#L70-L76)):**

```cpp
// Non-resonant (standard OP)
fully_integrate_d2Phi_a_domega_drho_in_rho(
  test_ergs, s, &SolarModel::Gamma_TP, 
  output_path + "TP.dat"
);

// Non-resonant (Rosseland)
fully_integrate_d2Phi_a_domega_drho_in_rho(
  test_ergs, s, &SolarModel::Gamma_TP_Rosseland, 
  output_path + "TP_Rosseland.dat"
);
```

**Integration Sequence:**
1. Radial integration over Sun's radius: $r \in [0, R_\odot]$
2. Energy integration per photon: $\omega \in [\omega_p(r), \omega_{\max}]$
3. Double integral: $\Phi = \int_{\omega_{\min}}^{\omega_{\max}} d\omega \int_0^{R_\odot} dr \cdot \frac{d^2\Phi}{d\omega dr}$

---

## 4. DATA COMPARISON & NORMALIZATION ISSUES

### 4.1 Flux Normalization Discrepancy

**From [plot_results.py](SolarAxionFlux/results/plot_results.py#L121):**

```python
# Hoof (OP-based) vs. Giannotti reference data
ax.plot(res6[:,0], res6[:,1]/1.0e10, 'k-', 
        label=r'TP (AGSS09)')

# Giannotti benchmark with correction factors
ax.plot(ref5[:,0], ref5[:,1]*4.0*1.4995, '-', 
        color='green', label=r'TP (Giannotti)')
```

**Correction Factor Breakdown:**
- Base Giannotti flux: $\Phi_{\text{Giannotti}}$ (dimensionless per $g_{a\gamma}^2$)
- Hoof OP flux: $\Phi_{\text{Hoof,OP}}/10^{10}$ (after coupling normalization)
- Effective correction: $4.0 \times 1.4995 = 5.998 \approx 6$

**Interpretation:**

Factor of $4.0$:
- Likely related to different definitions of $g_{a\gamma}$ coupling constant
- Possible: $(g_{a\gamma}^{\text{Giannotti}})^2 / (g_{a\gamma}^{\text{Hoof}})^2 = 4$
- Or magnetic field normalization difference

Factor of $1.4995 \approx 1.5 = 3/2$:
- Consistent with B-field averaging difference
- Giannotti: $\langle B_T^2 \rangle$ (unspecified averaging)
- Hoof: $\langle B^2 \rangle/3$ (isotropic)
- Ratio: $\frac{1}{3} / \frac{1}{3} \cdot f_{B} = 1.5$ suggests Giannotti uses different field profile

### 4.2 Energy Range Differences

**Typical Integration Ranges:**

| Implementation | $\omega_{\min}$ | $\omega_{\max}$ | Notes |
|---------------|--------------|--------------|-------|
| Hoof (OP) | $\omega_p(r)$ (energy-dependent) | ~19 keV | Full spectrum |
| Hoof (Rosseland) | $\omega_p(r)$ (energy-dependent) | ~19 keV | Full spectrum |
| Giannotti | Not specified | Likely up to 20 keV | Probably similar |

**Solar Interior Plasma Frequency:**
- At core ($r \approx 0$): $\omega_p \approx 0.5$ keV
- At photosphere: $\omega_p \ll 0.01$ keV
- Integration always includes core region where $\omega_p$ is significant

### 4.3 Data File Outputs

**Hoof Generates (from [tests.hpp](include/solaxflux/tests.hpp#L70-L111)):**
1. `TP.dat` - Standard OP-based spectrum
2. `TP_Rosseland.dat` - Rosseland-based spectrum
3. `TP_resonant_m131.dat`, `TP_resonant_m11.dat`, `TP_resonant_m0.dat` - Massive/resonant modes

**Giannotti Reference Data:**
- `2020_giannotti_TP.dat` - Reference flux (pre-computed)
- Used for validation and comparison

---

## 5. RESONANT (MASSIVE) TP CONTRIBUTIONS

### 5.1 Extension to Non-Zero Axion Masses

The Hoof code also implements **resonant TP production** for specific axion masses:

```cpp
double SolarModel::Gamma_TP_resonant_m131(double omega, double r);
double SolarModel::Gamma_TP_resonant_m11(double omega, double r);
double SolarModel::Gamma_TP_resonant_m0(double omega, double r);
```

**Resonance Structure:**
When $m_a \neq 0$, the denominator becomes:
$$(\Delta_p - \Delta_a)^2 + (\gamma/2)^2$$

Where $\Delta_a = -m_a^2/(2\omega)$ causes a **Breit-Wigner resonance** at:
$$\omega_{\text{res}}^2 \approx \omega_p^2 + m_a^2$$

**Specific Masses in Hoof:**
- $m_a = 0$ eV: Off-resonance (Giannotti case)
- $m_a = 11$ eV: Sub-keV resonance window
- $m_a = 131$ eV: Higher-energy resonance

### 5.2 Why Giannotti Focuses on $m_a = 0$

Giannotti et al. specifically study the **non-resonant zero-mass flux** because:
1. **Model-independent**: Valid for all ALP models with $m_a < 100$ eV
2. **Dominant at low energies**: Most photons below 1 keV couple off-resonance
3. **Experimental relevance**: IAXO and other helioscopes most sensitive in $\omega \sim 1-10$ keV

---

## 6. SUMMARY TABLE: GIANNOTTI VS. HOOF

| Parameter | Giannotti | Hoof (OP) | Hoof (Rosseland) |
|-----------|-----------|-----------|------------------|
| **Opacity Model** | Unspecified (likely OP) | Opacity Project (element-by-element) | Rosseland mean (energy-weighted) |
| **Opacity Energy Dependence** | Yes | Yes | No (radius only) |
| **Momentum Transfer Formula** | Simplified ($\Delta_p^2 = \omega_p^4/4\omega^2$) | Full ($\Delta_P^2 = \omega^2[\sqrt{1-\omega_p^2/\omega^2}-1]^2$) | Same as OP |
| **Magnetic Field Averaging** | Unspecified | Isotropic ($B^2/3$) | Same as OP |
| **Plasma Frequency Cutoff** | Soft (implicit) | Hard ($\omega > \omega_p$ enforced) | Same as OP |
| **Numeric Stability** | Standard | High (expm1, pow_2, etc.) | Same as OP |
| **Computational Cost** | N/A (reference data) | High (energy-dependent opacity lookups) | Low (single Rosseland value per radius) |
| **Flux Normalization** | Baseline (1.0) | 10^{-10} × (Hoof result) | Same as OP |
| **Correction Factor to Match** | — | × 4.0 × 1.4995 ≈ 6 | Same ratio |
| **Applicable Mass Range** | $m_a \lesssim 100$ eV (non-resonant) | $m_a = 0$ (also resonant modes available) | $m_a = 0$ |

---

## 7. PHYSICAL ISSUES & RECOMMENDATIONS

### 7.1 Outstanding Physics Questions

1. **Opacity Model Dependence**:
   - How sensitive is TP flux to opacity code choice (OP vs. LEDCOP vs. ATOMIC)?
   - Giannotti results may use different opacity treatment
   - **Recommendation**: Cross-check with LEDCOP/ATOMIC codes

2. **Momentum Transfer Formula**:
   - Hoof's full formula reduces to Giannotti's simplified form in certain limits
   - Low-energy discrepancy: Where exactly does $\sqrt{1-x}$ expansion break down?
   - **Recommendation**: Verify low-energy ($\omega < 1$ keV) behavior analytically

3. **Magnetic Field Profile**:
   - Hoof's stratified model (radial + tachocline + outer) vs. Giannotti's unspecified field
   - Possible explanation for ~1.5 factor difference
   - **Recommendation**: Test alternative B-field models in Hoof code

4. **Thermal Damping Rate**:
   - Both use $\gamma = -\text{expm1}(-\omega/T) \cdot \kappa$
   - But definitions of $\kappa$ differ (frequency-dependent vs. mean)
   - At what energies does this dominate the difference?

### 7.2 Numerical Issues to Monitor

**Hoof Code Robustness:**
1. ✓ Proper use of `expm1` avoids catastrophic cancellation
2. ✓ Early plasma frequency exit prevents NaN
3. ✓ Isotropic field averaging mathematically well-defined
4. ⚠️ Momentum transfer formula valid only for $\omega > \omega_p$; check boundary

**Potential Problems:**
- Very low energies ($\omega \to 0$): Behavior of $\sqrt{1-\omega_p^2/\omega^2}$
- Very high energies ($\omega \gg 20$ keV): Integration cutoff effects
- Interfaces between solar zones: Discontinuities in B-field

---

## 8. COMPARATIVE ADVANTAGES

### Giannotti Paper Advantages:
1. **Theoretical clarity**: Explicit focus on non-resonant physics
2. **Model independence**: No opacity code specification (applicable broadly)
3. **Coherent picture**: Large-scale magnetic field conversion mechanism

### Hoof Code Advantages:
1. **Implementation quality**: High numerical stability, well-tested
2. **Flexibility**: Multiple opacity codes (OP, LEDCOP, ATOMIC, OPAS)
3. **Extensibility**: Includes resonant (massive) contributions beyond Giannotti
4. **Practical tool**: Easily compute spectra for arbitrary ALPs

---

## 9. RECONCILIATION STRATEGY

To fully reconcile Giannotti vs. Hoof:

1. **Step 1**: Verify opacity model
   - Extract/compute Rosseland opacity used by Giannotti
   - Recompute Hoof TP with that exact opacity

2. **Step 2**: Check magnetic field assumptions
   - Determine Giannotti's B-field profile from paper
   - Implement in Hoof and compare factor

3. **Step 3**: Compare low-energy behavior
   - Extract limiting form of momentum transfer in both
   - Verify Taylor expansion matches

4. **Step 4**: Full spectral comparison
   - Plot Giannotti vs. Hoof (both OP and Rosseland)
   - Identify energy ranges with largest discrepancies

---

## References

- **Giannotti et al.** (2020): "Production of axionlike particles from photon conversions in large-scale solar magnetic fields", *Phys. Rev. D* **102**, 123024
- **Hoof et al.** (2021): "Quantifying uncertainties in the solar axion flux and their impact on determining axion model parameters", *JCAP* **09**, 006
- **Carenza et al.** (2022): "Axion-electron scattering in astrophysical environments", *arXiv:2201.07803*

---

## Appendix: Formula Cross-Reference

### Transverse Polarization Production Rate

**Giannotti Form (Equation in Paper):**
$$\dot{n}_\gamma \to a = N_{\text{pol}} \cdot |\text{Ampl}|^2 \cdot \rho(\omega) \cdot [1 + n_B(\omega)]$$

**Where amplitude squared is proportional to:**
$$\frac{\Delta_T^2 \cdot \gamma}{[\Delta_p - \Delta_a]^2 + [\gamma/2]^2}$$

**Hoof Implementation (Dimensionless):**
$$\Gamma_{\text{TP}} = 2 \cdot \frac{\gamma \cdot \Delta_T^2}{[\Delta_P^2 + (\gamma/2)^2] \cdot (e^{\omega/T}-1)}$$

The factor $(e^{\omega/T}-1)^{-1} = 1/(e^{\omega/T}-1)$ appears because Hoof works with production rate per unit volume per unit energy, requiring division by Bose-Einstein factor.

