//! Adsorption profiles and isotherms.
use super::functional::{HelmholtzEnergyFunctional, DFT};
use super::profile::DFTSpecifications;
use super::solver::DFTSolver;
use feos_core::{
    Contributions, EosError, EosResult, EosUnit, EquationOfState, StateBuilder, VLEOptions,
};
use ndarray::{arr1, Array1, Dimension, Ix1, Ix3};
use quantity::{QuantityArray1, QuantityArray2, QuantityScalar};
use std::rc::Rc;

mod external_potential;
mod fea_potential;
mod pore;
pub use external_potential::{ExternalPotential, FluidParameters};
pub use pore::{Pore1D, Pore3D, PoreProfile, PoreProfile1D, PoreProfile3D, PoreSpecification};

const MAX_ITER_ADSORPTION_EQUILIBRIUM: usize = 50;
const TOL_ADSORPTION_EQUILIBRIUM: f64 = 1e-8;

/// Possible inputs for the quantity grid of adsorption isotherms/isosteres.
pub enum QuantitySpecification<U> {
    /// Specify the minimal and maximal pressure, and the number of points
    Qlim {
        q_min: QuantityScalar<U>,
        q_max: QuantityScalar<U>,
        points: usize,
    },
    /// Provide a custom array of pressure points.
    Qvec(QuantityArray1<U>),
}

impl<U: EosUnit> QuantitySpecification<U>
where
    QuantityScalar<U>: std::fmt::Display,
{
    fn q_min_max(&self) -> (QuantityScalar<U>, QuantityScalar<U>) {
        match self {
            Self::Qlim {
                q_min,
                q_max,
                points: _,
            } => (*q_min, *q_max),
            Self::Qvec(quantity) => (quantity.get(0), quantity.get(quantity.len() - 1)),
        }
    }

    fn to_vec(&self) -> EosResult<QuantityArray1<U>> {
        Ok(match self {
            Self::Qlim {
                q_min,
                q_max,
                points,
            } => QuantityArray1::linspace(*q_min, *q_max, *points)?,
            Self::Qvec(quantity) => quantity.clone(),
        })
    }

    fn equilibrium<D: Dimension, F: HelmholtzEnergyFunctional>(
        &self,
        equilibrium: &Adsorption<U, D, F>,
    ) -> EosResult<(QuantityArray1<U>, QuantityArray1<U>)>
    where
        D::Larger: Dimension<Smaller = D>,
    {
        let p_eq = equilibrium.pressure().get(0);
        match self {
            Self::Qlim {
                q_min,
                q_max,
                points,
            } => Ok((
                QuantityArray1::linspace(*q_min, p_eq, points / 2)?,
                QuantityArray1::linspace(*q_max, p_eq, points - points / 2)?,
            )),
            Self::Qvec(pressure) => {
                let index = (0..pressure.len()).find(|&i| pressure.get(i) > p_eq);
                Ok(if let Some(index) = index {
                    (
                        QuantityArray1::from_shape_fn(index + 1, |i| {
                            if i == index {
                                p_eq
                            } else {
                                pressure.get(i)
                            }
                        }),
                        QuantityArray1::from_shape_fn(pressure.len() - index + 1, |i| {
                            if i == pressure.len() - index {
                                p_eq
                            } else {
                                pressure.get(pressure.len() - i - 1)
                            }
                        }),
                    )
                } else {
                    (pressure.clone(), arr1(&[]) * U::reference_pressure())
                })
            }
        }
    }
}

/// Container structure for the calculation of adsorption isotherms.
pub struct Adsorption<U, D: Dimension, F>(pub Vec<EosResult<PoreProfile<U, D, F>>>, usize);

/// Container structure for adsorption isotherms in 1D pores.
pub type Adsorption1D<U, F> = Adsorption<U, Ix1, F>;
/// Container structure for adsorption isotherms in 3D pores.
pub type Adsorption3D<U, F> = Adsorption<U, Ix3, F>;

impl<U: EosUnit, D: Dimension, F: HelmholtzEnergyFunctional> Adsorption<U, D, F>
where
    QuantityScalar<U>: std::fmt::Display,
    D::Larger: Dimension<Smaller = D>,
{
    /// Calculate an adsorption isotherm (starting at low pressure)
    pub fn adsorption_isotherm<S: PoreSpecification<U, D, F>>(
        functional: &Rc<DFT<F>>,
        temperature: QuantityScalar<U>,
        pressure: &QuantitySpecification<U>,
        pore: &S,
        molefracs: Option<&Array1<f64>>,
        solver: Option<&DFTSolver>,
    ) -> EosResult<Adsorption<U, D, F>> {
        let pressure = pressure.to_vec()?;
        Self::isotherm(functional, temperature, &pressure, pore, molefracs, solver)
    }

    /// Calculate an desorption isotherm (starting at high pressure)
    pub fn desorption_isotherm<S: PoreSpecification<U, D, F>>(
        functional: &Rc<DFT<F>>,
        temperature: QuantityScalar<U>,
        pressure: &QuantitySpecification<U>,
        pore: &S,
        molefracs: Option<&Array1<f64>>,
        solver: Option<&DFTSolver>,
    ) -> EosResult<Adsorption<U, D, F>> {
        let pressure = pressure.to_vec()?;
        let pressure =
            QuantityArray1::from_shape_fn(pressure.len(), |i| pressure.get(pressure.len() - i - 1));
        let isotherm = Self::isotherm(functional, temperature, &pressure, pore, molefracs, solver)?;
        Ok(Adsorption(
            isotherm.0.into_iter().rev().collect(),
            functional.components(),
        ))
    }

    /// Calculate an equilibrium isotherm
    pub fn equilibrium_isotherm<S: PoreSpecification<U, D, F>>(
        functional: &Rc<DFT<F>>,
        temperature: QuantityScalar<U>,
        pressure: &QuantitySpecification<U>,
        pore: &S,
        molefracs: Option<&Array1<f64>>,
        solver: Option<&DFTSolver>,
    ) -> EosResult<Adsorption<U, D, F>> {
        let (p_min, p_max) = pressure.q_min_max();
        let equilibrium = Self::phase_equilibrium(
            functional,
            temperature,
            p_min,
            p_max,
            pore,
            molefracs,
            solver,
            VLEOptions::default(),
        );
        if let Ok(equilibrium) = equilibrium {
            let pressure = pressure.equilibrium(&equilibrium)?;
            let adsorption = Self::isotherm(
                functional,
                temperature,
                &pressure.0,
                pore,
                molefracs,
                solver,
            )?
            .0;
            let desorption = Self::isotherm(
                functional,
                temperature,
                &pressure.1,
                pore,
                molefracs,
                solver,
            )?
            .0;
            Ok(Adsorption(
                adsorption
                    .into_iter()
                    .chain(desorption.into_iter().rev())
                    .collect(),
                functional.components(),
            ))
        } else {
            let adsorption = Self::adsorption_isotherm(
                functional,
                temperature,
                pressure,
                pore,
                molefracs,
                solver,
            )?;
            let desorption = Self::desorption_isotherm(
                functional,
                temperature,
                pressure,
                pore,
                molefracs,
                solver,
            )?;
            let omega_a = adsorption.grand_potential();
            let omega_d = desorption.grand_potential();
            let is_ads = Array1::from_shape_fn(adsorption.0.len(), |i| {
                omega_d.get(i).is_nan() || omega_a.get(i) < omega_d.get(i)
            });
            let profiles = is_ads
                .into_iter()
                .zip(adsorption.0.into_iter())
                .zip(desorption.0.into_iter())
                .map(|((is_ads, a), d)| if is_ads { a } else { d })
                .collect();
            Ok(Adsorption(profiles, functional.components()))
        }
    }

    fn isotherm<S: PoreSpecification<U, D, F>>(
        functional: &Rc<DFT<F>>,
        temperature: QuantityScalar<U>,
        pressure: &QuantityArray1<U>,
        pore: &S,
        molefracs: Option<&Array1<f64>>,
        solver: Option<&DFTSolver>,
    ) -> EosResult<Adsorption<U, D, F>> {
        let moles =
            functional.validate_moles(molefracs.map(|x| x * U::reference_moles()).as_ref())?;
        let mut profiles: Vec<EosResult<PoreProfile<U, D, F>>> = Vec::with_capacity(pressure.len());

        // Calculate the external potential once
        let mut bulk = StateBuilder::new(functional)
            .temperature(temperature)
            .pressure(pressure.get(0))
            .moles(&moles)
            .build()?;
        if functional.components() > 1 && !bulk.is_stable(VLEOptions::default())? {
            bulk = bulk
                .tp_flash(None, VLEOptions::default(), None)?
                .vapor()
                .clone();
        }
        let external_potential = pore
            .initialize(&bulk, None, None)?
            .profile
            .external_potential;

        for i in 0..pressure.len() {
            let mut bulk = StateBuilder::new(functional)
                .temperature(temperature)
                .pressure(pressure.get(i))
                .moles(&moles)
                .build()?;
            if functional.components() > 1 && !bulk.is_stable(VLEOptions::default())? {
                bulk = bulk
                    .tp_flash(None, VLEOptions::default(), None)?
                    .vapor()
                    .clone();
            }
            let mut p = pore.initialize(&bulk, Some(&external_potential), None)?;
            let p2 = p.clone();
            if let Some(Ok(l)) = profiles.last() {
                p.profile.density = l.profile.density.clone();
            }
            profiles.push(p.solve(solver).or_else(|_| p2.solve(solver)));
        }

        Ok(Adsorption(profiles, functional.components()))
    }

    /// Calculate the phase transition from an empty to a filled pore.
    pub fn phase_equilibrium<S: PoreSpecification<U, D, F>>(
        functional: &Rc<DFT<F>>,
        temperature: QuantityScalar<U>,
        p_min: QuantityScalar<U>,
        p_max: QuantityScalar<U>,
        pore: &S,
        molefracs: Option<&Array1<f64>>,
        solver: Option<&DFTSolver>,
        options: VLEOptions,
    ) -> EosResult<Adsorption<U, D, F>> {
        let moles =
            functional.validate_moles(molefracs.map(|x| x * U::reference_moles()).as_ref())?;

        // calculate density profiles for the minimum and maximum pressure
        let vapor_bulk = StateBuilder::new(functional)
            .temperature(temperature)
            .pressure(p_min)
            .moles(&moles)
            .vapor()
            .build()?;
        let liquid_bulk = StateBuilder::new(functional)
            .temperature(temperature)
            .pressure(p_max)
            .moles(&moles)
            .liquid()
            .build()?;

        let mut vapor = pore.initialize(&vapor_bulk, None, None)?.solve(None)?;
        let mut liquid = pore.initialize(&liquid_bulk, None, None)?.solve(solver)?;

        // calculate initial value for the molar gibbs energy
        let nv = vapor.profile.bulk.density
            * (vapor.profile.moles() * vapor.profile.bulk.molar_volume(Contributions::Total)).sum();
        let nl = liquid.profile.bulk.density
            * (liquid.profile.moles() * liquid.profile.bulk.molar_volume(Contributions::Total))
                .sum();
        let f = |s: &mut PoreProfile<U, D, F>, n: QuantityScalar<U>| -> EosResult<_> {
            Ok(s.grand_potential.unwrap()
                + s.profile.bulk.molar_gibbs_energy(Contributions::Total) * n)
        };
        let mut g = (f(&mut liquid, nl)? - f(&mut vapor, nv)?) / (nl - nv);

        // update filled pore with limited step size
        let mut bulk = StateBuilder::new(functional)
            .temperature(temperature)
            .pressure(p_max)
            .moles(&moles)
            .vapor()
            .build()?;
        let g_liquid = liquid.profile.bulk.molar_gibbs_energy(Contributions::Total);
        let steps = (10.0 * (g - g_liquid)).to_reduced(g_liquid)?.abs().ceil() as usize;
        let delta_g = (g - g_liquid) / steps as f64;
        for i in 1..=steps {
            let g_i = g_liquid + i as f64 * delta_g;
            bulk = bulk.update_gibbs_energy(g_i)?;
            liquid = liquid.update_bulk(&bulk).solve(solver)?;
        }

        for _ in 0..options.max_iter.unwrap_or(MAX_ITER_ADSORPTION_EQUILIBRIUM) {
            // update empty pore
            vapor = vapor.update_bulk(&bulk).solve(None)?;

            // update filled pore
            liquid = liquid.update_bulk(&bulk).solve(solver)?;

            // calculate moles
            let nv = vapor.profile.bulk.density
                * (vapor.profile.moles() * vapor.profile.bulk.molar_volume(Contributions::Total))
                    .sum();
            let nl = liquid.profile.bulk.density
                * (liquid.profile.moles() * liquid.profile.bulk.molar_volume(Contributions::Total))
                    .sum();

            // check for a trivial solution
            if nl.to_reduced(nv)? - 1.0 < 1e-5 {
                return Err(EosError::TrivialSolution);
            }

            // Newton step
            let delta_g =
                (vapor.grand_potential.unwrap() - liquid.grand_potential.unwrap()) / (nv - nl);
            if delta_g.to_reduced(U::reference_molar_energy())?.abs()
                < options.tol.unwrap_or(TOL_ADSORPTION_EQUILIBRIUM)
            {
                return Ok(Adsorption(
                    vec![Ok(vapor), Ok(liquid)],
                    functional.components(),
                ));
            }
            g += delta_g;

            // update bulk phase
            bulk = bulk.update_gibbs_energy(g)?;
        }
        Err(EosError::NotConverged(
            "Adsorption::phase_equilibrium".into(),
        ))
    }

    pub fn isostere<S: PoreSpecification<U, D, F>>(
        &self,
        functional: &Rc<DFT<F>>,
        initial_pressure: QuantityScalar<U>,
        temperature: &QuantitySpecification<U>,
        pore: &S,
        molefracs: Option<&Array1<f64>>,
        solver: Option<&DFTSolver>,
    ) -> EosResult<Adsorption<U, D, F>> {
        let temperature = temperature.to_vec()?;
        let moles =
            functional.validate_moles(molefracs.map(|x| x * U::reference_moles()).as_ref())?;
        let mut profiles: Vec<EosResult<PoreProfile<U, D, F>>> =
            Vec::with_capacity(temperature.len());

        // Create the initial bulk state at the start of the isostere
        let mut initial_bulk = StateBuilder::new(functional)
            .temperature(temperature.get(0))
            .pressure(initial_pressure)
            .moles(&moles)
            .build()?;
        if functional.components() > 1 && !initial_bulk.is_stable(VLEOptions::default())? {
            initial_bulk = initial_bulk
                .tp_flash(None, VLEOptions::default(), None)?
                .vapor()
                .clone();
        }

        // Initial PoreProfile for given mu,V,T
        let initial_profile = pore.initialize(&initial_bulk, None, None)?.solve(solver)?;

        let specification = DFTSpecifications::moles_from_profile(&initial_profile.profile)?;

        profiles.push(Ok(initial_profile.clone()));

        let partial_density = initial_profile.profile.moles() / initial_profile.profile.volume();

        for i in 1..temperature.len() {
            let dummy_state = StateBuilder::new(functional)
                .temperature(temperature.get(i))
                .partial_density(&partial_density)
                .build()?;

            let mut p = pore
                .initialize(&dummy_state, None, Some(specification.clone()))?
                .solve(solver)?;
            p.profile
                .bulk
                .update_chemical_potential(&p.profile.chemical_potential)?;

            profiles.push(Ok(p));
        }

        Ok(Adsorption(profiles, functional.components()))
    }

    pub fn pressure(&self) -> QuantityArray1<U> {
        QuantityArray1::from_shape_fn(self.0.len(), |i| match &self.0[i] {
            Ok(p) => {
                if p.profile.bulk.eos.components() > 1
                    && !p.profile.bulk.is_stable(VLEOptions::default()).unwrap()
                {
                    p.profile
                        .bulk
                        .tp_flash(None, VLEOptions::default(), None)
                        .unwrap()
                        .vapor()
                        .pressure(Contributions::Total)
                } else {
                    p.profile.bulk.pressure(Contributions::Total)
                }
            }
            Err(_) => f64::NAN * U::reference_pressure(),
        })
    }

    pub fn molar_gibbs_energy(&self) -> QuantityArray1<U> {
        QuantityArray1::from_shape_fn(self.0.len(), |i| match &self.0[i] {
            Ok(p) => {
                if p.profile.bulk.eos.components() > 1
                    && !p.profile.bulk.is_stable(VLEOptions::default()).unwrap()
                {
                    p.profile
                        .bulk
                        .tp_flash(None, VLEOptions::default(), None)
                        .unwrap()
                        .vapor()
                        .molar_gibbs_energy(Contributions::Total)
                } else {
                    p.profile.bulk.molar_gibbs_energy(Contributions::Total)
                }
            }
            Err(_) => f64::NAN * U::reference_molar_energy(),
        })
    }

    pub fn adsorption(&self) -> QuantityArray2<U> {
        QuantityArray2::from_shape_fn((self.1, self.0.len()), |(j, i)| match &self.0[i] {
            Ok(p) => p.profile.moles().get(j),
            Err(_) => f64::NAN * U::reference_moles() / U::reference_length().powi(2),
        })
    }

    pub fn total_adsorption(&self) -> QuantityArray1<U> {
        QuantityArray1::from_shape_fn(self.0.len(), |i| match &self.0[i] {
            Ok(p) => p.profile.total_moles(),
            Err(_) => f64::NAN * U::reference_moles() / U::reference_length().powi(2),
        })
    }

    pub fn grand_potential(&self) -> QuantityArray1<U> {
        QuantityArray1::from_shape_fn(self.0.len(), |i| match &self.0[i] {
            Ok(p) => p.grand_potential.unwrap(),
            Err(_) => f64::NAN * U::reference_surface_tension(),
        })
    }
}
