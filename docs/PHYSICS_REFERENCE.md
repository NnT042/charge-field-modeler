# Charge Field Theory — Physics Reference for CFM Development

This document summarizes the key physics concepts from Miles Mathis's charge field theory that are directly relevant to implementing the Charge Field Modeler. It is not a comprehensive overview of the theory — only the parts that translate into simulation rules, data structures, and testable predictions.

Developers working on this project should read this file before writing any simulation code.

---

## 1. Foundational Principles

### No Point Particles
Every particle has a real, physical radius. There are no point particles and no zero-dimensional objects. The base photon (B-photon) has a calculable radius approximately equal to G × proton radius ≈ 6.67 × 10⁻²⁶ m.

### No Attraction
All charge is repulsive. There is no attractive electromagnetic force. What appears as attraction is always the result of charge field geometry — a particle being pushed toward another by the surrounding field, not pulled by the other particle. Gravity at the quantum level is real and is a function of radius (treated mathematically as expansion of the sphere).

### Charge = Real Photon Emission
The charge field is composed of real, physical photons (B-photons) traveling at c. These are not virtual particles. Charge is mediated by straight bombardment — photons hitting things and transferring momentum. The strength of a particle's charge is a function of its surface area (how many B-photons it can emit per unit time).

---

## 2. Stacked Spins

### The Four Spin Types

A spherical particle can carry up to four independent, simultaneous spins:

1. **Axial (a-spin):** Rotation about the particle's own axis. The radius of this spin equals the particle's physical radius r.

2. **X-spin:** End-over-end tumble in the x-plane. Amplitude = 2r. This spin is outside the gyroscopic influence of the axial spin because it is exactly double the radius — the minimum distance needed to avoid interference.

3. **Y-spin:** End-over-end in the y-plane. Amplitude = 4r. Outside the gyroscopic influence of the x-spin for the same doubling reason.

4. **Z-spin:** End-over-end in the z-plane. Amplitude = 8r. Outside the gyroscopic influence of the y-spin.

### Why Doubling (Quantization)

The doubling at each level is the mechanical cause of quantization. Each spin level must be exactly double the previous to avoid gyroscopic interference. Fractional values are impossible — you cannot have a 1.5r spin because it would interfere with the axial spin. This is why energy is quantized: E = hν, and ν is quantized because the wavelength is quantized, and the wavelength is quantized because it is a stacked spin with mandatory doubling.

### Wavelength from Spin

When a particle with stacked spins has linear motion, the outermost spin traces out the measurable wavelength. For a particle with all four spins and radius r, the z-spin has amplitude 8r, so the electrical wavelength = 8r. To find the particle's actual radius from a measured wavelength: divide by 8.

### Chirality

Each spin level can be clockwise (+) or counterclockwise (-). Per Mathis (quark.html): "I assign clockwise to +a, and I do this simply to keep the right-hand rule." Our implementation matches: `DQuat::from_axis_angle` with positive angle produces CW rotation when viewed from behind the direction of travel.

The combination of chiralities across all levels determines the particle type (proton vs neutron vs anti-proton vs anti-neutron) at the baryon scale.

### Angular Velocity Convention

Our sliders use `v` in natural units where `|v| = 1.0` is c (speed of light). The simulation computes `ω = v / orbit_radius`, so outer levels rotate slower. This uses the standard v=ωr formula. Mathis corrects this to v=ω/r for quantum particles, but the correction only affects absolute SI rates (absorbed by the time_scale slider) — the relative scaling between levels is preserved either way.

### Pi = 4 (Kinematic vs Geometric)

Mathis argues kinematic circumference is 8r, not 2πr. We use both correctly:
- **Energy formula**: uses 8r per level (kinematic — measures energy from motion through the orbit)
- **Wavelength readout**: uses 2πr (geometric — the spatial period of the traced wave)
- **advance()**: uses standard angular rate (the ratio between levels is unchanged by the π correction)

---

## 3. Particle Hierarchy

### Photon Tier (Levels 1-4)
A base photon with 1-4 spin levels. Travels at c. The z-spin creates the electrical wavelength. Higher-energy photons have more spin levels active. A 4-level photon is a high-energy photon (UV range and above, depending on base radius).

### Electron Tier (Levels 5-8)
The 4-level photon structure is treated as a new effective sphere, and a second set of four spins is stacked on top. Level 8 (z₂) completes the "non-spinning electron" — all photon-tier spins saturated, but no baryon-level spin yet. From Mathis (mit.pdf): spin levels exist structurally but can be "torpid" (not energized), which maps to our ω=0 mechanic.

### Baryon Tier (Levels 9-12)
Level 9 is the **electron at rest** — it adds baryon-tier axial spin (energy = 9 in natural units). From Mathis (photon.html): "The electron is spinning axially; the proton is spinning axially plus x, y, and z." The proton adds x, y, z spins at the baryon tier (levels 10-12), giving energy 16385. Ratio: 16385/9 = 1820.56 ≈ the measured nucleon/electron mass ratio.

The difference between proton and neutron is the combination of spin chiralities:

**Proton spin states (from QCD overhaul paper):**
- +a+x+y+z
- -a-x-y-z
- -a+x+y+z
- +a-x-y-z

These combinations allow B-photon emission to escape through the stacked spins. The emission exits in the same state it was emitted from the particle surface — same direction, same chirality. The proton is therefore charged (emitting).

**Neutron spin states:**
- -a-x-y+z
- +a+x+y-z
- -a+x+y-z
- +a-x-y+z

These combinations trap B-photon emission within the stacked spins. The emission cannot escape and returns with reversed spin, canceling its own energy. The neutron is therefore neutral (not emitting).

**Anti-proton spin states:**
- +a-x+y-z
- -a+x-y+z
- +a+x-y+z
- -a-x+y-z

Emission escapes but upside-down (reversed chirality relative to proton emission).

**Anti-neutron spin states:**
- +a-x+y+z
- +a+x-y-z
- -a-x+y+z
- -a+x-y-z

Emission is trapped but does not fully cancel — spin energy of the trapped photons is not reversed, leading to a slight mass difference from the neutron.

### Uberon Tier (Levels 13-15)
Levels above the baryon z-spin. From Mathis (meson.html): "'Meson' is not a logical term for a particle above the baryon mass, so I will call these particles 'uberons.'" Level 13 is the proton engine starter (baryon-tier axial precession). Level 14 is D meson territory (~1860 MeV). Level 15 is unstable uberon. These are collision-only transients — "The proton cannot add extra spins above the z-spin without creating instability." Level 16 (z₄, amplitude 4096) is excluded from the simulation: it overflows the display and has no stable physical analogue.

---

## 3b. Energy from Stacked Spins (The 1820 Derivation)

Source: "Unifying the Electron and Proton" (Mathis, 2008)

The electron and proton are the same fundamental particle. The difference is spin: the electron at rest has only axial spin, while the proton has all four spins at the baryon tier. The mass ratio between them (1820) is derived purely from the energy added by each spin level.

### The Formula

Start with a non-spinning particle at energy 1 (rest/linear energy). Add axial spin — using tangential velocity with π=4 (Mathis's kinematic correction), the circumference is 8r, giving spin energy 8. Total with axial spin: **[1 + 8] = 9**.

For each subsequent spin level, the energy multiplies because the spin is orthogonal, but is divided by powers of 2 because each end-over-end spin is in the forward direction only half the time (and this compounds across orthogonal planes):

| Spin Level | Energy Term | Total Energy |
|------------|-------------|--------------|
| Axial only | [1 + 8] | 9 |
| + X-spin | [1 + (8 × 16)/2] | 65 |
| + Y-spin | [1 + (8 × 16 × 32)/2²] | 1,025 |
| + Z-spin | [1 + (8 × 16 × 32 × 64)/2⁴] | 16,385 |

Where the multiplied terms come from 8r at each level (r = 1, 2, 4, 8 → 8×1=8, 8×2=16, 8×4=32, 8×8=64).

### The Ratio

16,385 / 9 = **1820.56**

This matches the measured nucleon/electron mass ratio. The electron at rest is a proton stripped of its outer spins (x, y, z at the baryon tier), retaining only the axial spin.

### Implications for the Simulation

- **Geometric amplitudes** (1, 2, 4, 8, 16...) define the actual orbital paths — use these for position computation and path tracing.
- **Energy values** from the formula above define the mass/energy — use these for the readout panel display.
- The two are conceptually separate: the radius doubling tells you *where* the particle goes; the energy formula tells you *how much energy* that motion represents.
- A "moving electron" gains energy from its motion, acquiring an x-spin (energy ~65/9 = 7.2 times rest). This is what creates the wave characteristic of a moving electron — it literally gains a second spin from collisions with field photons.
- An electron with no z-spin is a meson. An electron with all stable spins is a baryon.

---

## 4. Charge Channeling (Proton as Faraday Disc)

The proton recycles the charge field like a Faraday disc motor:
- Charge (B-photons) enters at the poles
- Charge is channeled through the spinning structure
- Charge exits at the equator, compressed into a plane by the high spin rate

This is why Miles draws protons as discs in nuclear diagrams — the disc represents the charge emission profile, not the particle's shape. The particle is still a sphere, but its charge emission is planar.

The emission exits at approximately **30° north and south** of the equator:
- Photon charge at +30° angle
- Antiphoton charge at -30° angle
- These average at the equator for net emission

This 30° angle arises from the geometry of how emission navigates through the z-spin boundary. It is a prediction of the spin geometry and should emerge from a faithful simulation without being hard-coded.

The neutron cannot channel charge this way — the opposing spins block the emission and redirect it back out through the poles. This is the mechanical basis for charge neutrality.

---

## 5. The Ambient Charge Field

### Composition Near Earth
The ambient charge field near Earth's surface is approximately:
- **2/3 photons** (left-spinning / CW)
- **1/3 antiphotons** (right-spinning / CCW)

This asymmetry is caused by the Sun, which spins in one direction and tends to convert antiphotons to photons as it channels the charge field.

### Field Direction
In the vicinity of matter, protons set the net charge field direction. The total charge field sums in the plane of the protons' equators (x,y plane in standard orientation). The field cannot sum along the polar axis because protons emit along their equators.

### Implications for Simulation
The "Earth ground level" preset should use:
- 67% CW / 33% CCW photon ratio
- A net field direction in one plane (representing the summed proton emission direction of nearby matter)
- Thermal energy distribution

### Chirality-Dependent Collisions (Not Yet Implemented)

Per Mathis, charge photon collisions are mechanically spin-dependent. From neut2.pdf: "collisions affect only the outer spins... The spins offset, like little gears." Two distinct outcomes depending on chirality match:

- **Same chirality, opposing travel** (double-negative): spins augment. Energy transfers from linear to spin. Softer bounce, photon exits slower with more spin energy. (Inverse Compton scattering.)
- **Opposite chirality, opposing travel** (single-negative): spins cancel. Energy transfers from spin to linear. Harder bounce, photon exits faster with less spin. (Compton scattering.)

From halbach.pdf: "think of them as cogs — they don't match, they catch, and force one another apart... repulsion implies a real loss of local spin. Energy has transferred from spin energy to linear."

The exit angle is affected because collisions are tangential edge hits. Same-spin deflects along the spin direction; opposite-spin deflects against it.

Our field sim currently uses pure specular reflection (no chirality dependence). Implementing this requires: (1) reading the existing chirality bit during collision, (2) modifying exit velocity based on match/mismatch, (3) coupling collision statistics back into the spin stack to allow outer-level velocity to equilibrate.

---

## 6. Mass, Energy, and the Neutron-Proton Mass Difference

The neutron has more mass than the proton because it traps its B-photon emission internally. Since energy = mass (E = mc²), the trapped emission adds mass. The mass equivalence of a baryon's emission field is approximately 2.3 × 10⁻³⁰ kg (the neutron-proton mass difference).

The anti-neutron has a slightly different mass than the neutron because its trapped emission does not fully cancel. Clockwise photons meeting clockwise photons don't cancel spin — they double it. So the trapped field in an anti-neutron retains spin energy that either adds to or subtracts from the particle's total mass, depending on the specific spin combination.

---

## 7. Key Predictions to Test in Simulation

1. **Proton charge emission at ~30° from equator** — should emerge from spin geometry
2. **Neutron charge blocking** — field photons hitting a neutron should scatter without coherent equatorial emission
3. **Electron riding the charge current** — electrons placed in a directional charge field should be carried along, too large to pass through proton channels but pushed by the field
4. **Spontaneous H2 orientation** — two protons placed near each other in a charge field should orient with equators facing each other
5. **Alpha particle formation** — two protons + two neutrons should find a stable configuration matching the disc diagrams in the nuclear papers (neutrons between protons, charge channeled through interior holes)
6. **Gravity as field effect** — a density gradient in the charge field (more photons from one direction) should produce net motion of the focus particle toward lower density. This is a stretch goal for later phases.

---

## 8. References

All source material is in the project knowledge base. Key papers by document:
- **overhaul_of_QCD_by_Miles_Mathis.pdf** — Primary reference for stacked spins, particle types, spin state tables, beta decay
- **magmom.pdf** — Magnetic moment derivation, 30° emission angle, Earth field composition, proton/neutron charge recycling diagrams
- **nuclear.pdf** — Nuclear structure, disc diagrams, alpha particle configurations, element building
- **Electrical_Charge_by_Miles_Mathis.pdf** — Charge field fundamentals, permittivity as gravity, unit analysis
- **Weak_Interaction.pdf** — W/Z as stacked-spin baryons, meson equation, energy scaling
- **Unifying_the_Electron_and_the_Proton_by_Miles_Mathis.pdf** — Electron/proton relationship, unified field
- **drude.pdf** — Charge photons as real carriers in conductors, photon field mechanics
- **per4.pdf** — Periodic table element construction, charge channeling through nuclei, conduction vs magnetism
