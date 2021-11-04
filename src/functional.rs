use crate::convolver::Convolver;
use crate::functional_contribution::*;
use crate::ideal_chain_contribution::IdealChainContribution;
use crate::weight_functions::{WeightFunction, WeightFunctionInfo, WeightFunctionShape};
use feos_core::{
    Contributions, EosResult, EosUnit, EquationOfState, HelmholtzEnergy, HelmholtzEnergyDual,
    MolarWeight, StateHD,
};
use ndarray::*;
use num_dual::*;
use petgraph::graph::{Graph, UnGraph};
use petgraph::visit::EdgeRef;
use petgraph::Directed;
use quantity::{QuantityArray, QuantityArray1, QuantityScalar};
use std::ops::{AddAssign, MulAssign};
use std::rc::Rc;

/// Wrapper struct for the [HelmholtzEnergyFunctional] trait.
#[derive(Clone)]
pub struct DFT<T> {
    pub functional: T,
    pub component_index: Array1<usize>,
    pub m: Array1<f64>,
    pub ideal_chain_contribution: IdealChainContribution,
}

impl<T> DFT<T> {
    /// Create a new DFT struct for a homosegmented Helmholtz energy functional.
    pub fn new_homosegmented(functional: T, m: &Array1<f64>) -> Self {
        let component_index = Array1::from_shape_fn(m.len(), |i| i);
        Self::new(functional, &component_index, m)
    }

    /// Create a new DFT struct for a heterosegmented Helmholtz energy functional.
    pub fn new_heterosegmented(functional: T, component_index: &Array1<usize>) -> Self {
        let m = Array1::ones(component_index.len());
        Self::new(functional, component_index, &m)
    }

    /// Create a new DFT struct for a general Helmholtz energy functional.
    pub fn new(functional: T, component_index: &Array1<usize>, m: &Array1<f64>) -> Self {
        Self {
            functional,
            component_index: component_index.clone(),
            m: m.clone(),
            ideal_chain_contribution: IdealChainContribution::new(component_index, m),
        }
    }
}

impl<T: MolarWeight<U>, U: EosUnit> MolarWeight<U> for DFT<T> {
    fn molar_weight(&self) -> QuantityArray1<U> {
        self.functional.molar_weight()
    }
}

impl<T: HelmholtzEnergyFunctional> EquationOfState for DFT<T> {
    fn components(&self) -> usize {
        self.component_index[self.component_index.len() - 1] + 1
    }

    fn subset(&self, component_list: &[usize]) -> Self {
        self.functional.subset(component_list)
    }

    fn compute_max_density(&self, moles: &Array1<f64>) -> f64 {
        self.functional.compute_max_density(moles)
    }

    fn residual(&self) -> &[Box<dyn HelmholtzEnergy>] {
        unreachable!()
    }

    fn evaluate_residual<D: DualNum<f64>>(&self, state: &StateHD<D>) -> D
    where
        dyn HelmholtzEnergy: HelmholtzEnergyDual<D>,
    {
        self.functional
            .contributions()
            .iter()
            .map(|c| (c as &dyn HelmholtzEnergy).helmholtz_energy(state))
            .sum::<D>()
            + self.ideal_chain_contribution.helmholtz_energy(state)
    }

    fn evaluate_residual_contributions<D: DualNum<f64>>(
        &self,
        state: &StateHD<D>,
    ) -> Vec<(String, D)>
    where
        dyn HelmholtzEnergy: HelmholtzEnergyDual<D>,
    {
        let mut res: Vec<(String, D)> = self
            .functional
            .contributions()
            .iter()
            .map(|c| {
                (
                    c.to_string(),
                    (c as &dyn HelmholtzEnergy).helmholtz_energy(state),
                )
            })
            .collect();
        res.push((
            self.ideal_chain_contribution.to_string(),
            self.ideal_chain_contribution.helmholtz_energy(state),
        ));
        res
    }
}

/// A general Helmholtz energy functional.
pub trait HelmholtzEnergyFunctional: Sized {
    /// Return a slice of [FunctionalContribution]s.
    fn contributions(&self) -> &[Box<dyn FunctionalContribution>];

    /// Return a [DFT] for the specified subset of components.
    fn subset(&self, component_list: &[usize]) -> DFT<Self>;

    /// Return the maximum density in Angstrom^-3.
    ///
    /// This value is used as an estimate for a liquid phase for phase
    /// equilibria and other iterations. It is not explicitly meant to
    /// be a mathematical limit for the density (if those exist in the
    /// equation of state anyways).
    fn compute_max_density(&self, moles: &Array1<f64>) -> f64;

    /// Overwrite this, if the functional consists of heterosegmented chains.
    fn bond_lengths(&self, _temperature: f64) -> UnGraph<(), f64> {
        Graph::with_capacity(0, 0)
    }

    fn weight_functions(&self, temperature: f64) -> Vec<WeightFunctionInfo<f64>> {
        self.contributions()
            .iter()
            .map(|c| c.weight_functions(temperature))
            .collect()
    }
}

impl<T: HelmholtzEnergyFunctional> DFT<T> {
    pub fn grand_potential_density<U, D>(
        &self,
        temperature: QuantityScalar<U>,
        density: &QuantityArray<U, D::Larger>,
        convolver: &Rc<dyn Convolver<f64, D>>,
    ) -> EosResult<QuantityArray<U, D>>
    where
        U: EosUnit,
        D: Dimension,
        D::Larger: Dimension<Smaller = D>,
    {
        // Calculate residual Helmholtz energy density and functional derivative
        let t = temperature.to_reduced(U::reference_temperature()).unwrap();
        let rho = density.to_reduced(U::reference_density()).unwrap();
        let (mut f, dfdrho) = self.functional_derivative(t, &rho, convolver)?;

        // calculate the grand potential density
        for ((rho, dfdrho), &m) in rho.outer_iter().zip(dfdrho.outer_iter()).zip(self.m.iter()) {
            f -= &((&dfdrho + m) * &rho);
        }

        let bond_lengths = self.functional.bond_lengths(t);
        for segment in bond_lengths.node_indices() {
            let n = bond_lengths.neighbors(segment).count();
            f += &(&rho.index_axis(Axis(0), segment.index()) * (0.5 * n as f64));
        }

        Ok(f * t * U::reference_pressure())
    }

    fn intrinsic_helmholtz_energy_density<D, N>(
        &self,
        temperature: N,
        density: &Array<f64, D::Larger>,
        convolver: &Rc<dyn Convolver<N, D>>,
        contributions: Contributions,
    ) -> EosResult<Array<N, D>>
    where
        N: DualNum<f64> + ScalarOperand,
        dyn FunctionalContribution: FunctionalContributionDual<N>,
        D: Dimension,
        D::Larger: Dimension<Smaller = D>,
    {
        let density_dual = density.mapv(N::from);
        let weighted_densities = convolver.weighted_densities(&density_dual);
        let functional_contributions = self.functional.contributions();
        let mut helmholtz_energy_density: Array<N, D> = self
            .ideal_chain_contribution
            .calculate_helmholtz_energy_density(&density.mapv(N::from), contributions)?;
        for (c, wd) in functional_contributions.iter().zip(weighted_densities) {
            let nwd = wd.shape()[0];
            let ngrid = wd.len() / nwd;
            helmholtz_energy_density
                .view_mut()
                .into_shape(ngrid)
                .unwrap()
                .add_assign(&c.calculate_helmholtz_energy_density(
                    temperature,
                    wd.into_shape((nwd, ngrid)).unwrap().view(),
                )?);
        }
        Ok(helmholtz_energy_density * temperature)
    }

    pub fn entropy_density<D>(
        &self,
        temperature: f64,
        density: &Array<f64, D::Larger>,
        convolver: &Rc<dyn Convolver<Dual64, D>>,
        contributions: Contributions,
    ) -> EosResult<Array<f64, D>>
    where
        D: Dimension,
        D::Larger: Dimension<Smaller = D>,
    {
        let temperature_dual = Dual64::from(temperature).derive();
        let helmholtz_energy_density = self.intrinsic_helmholtz_energy_density(
            temperature_dual,
            density,
            convolver,
            contributions,
        )?;
        Ok(helmholtz_energy_density.mapv(|f| -f.eps[0]))
    }

    pub fn entropy_density_contributions<D>(
        &self,
        temperature: f64,
        density: &Array<f64, D::Larger>,
        convolver: &Rc<dyn Convolver<Dual64, D>>,
        contributions: Contributions,
    ) -> EosResult<Vec<Array<f64, D>>>
    where
        D: Dimension,
        D::Larger: Dimension<Smaller = D>,
        <D::Larger as Dimension>::Larger: Dimension<Smaller = D::Larger>,
    {
        let density_dual = density.mapv(Dual64::from);
        let temperature_dual = Dual64::from(temperature).derive();
        let weighted_densities = convolver.weighted_densities(&density_dual);
        let functional_contributions = self.functional.contributions();
        let mut helmholtz_energy_density: Vec<Array<Dual64, D>> =
            Vec::with_capacity(functional_contributions.len() + 1);
        helmholtz_energy_density.push(
            self.ideal_chain_contribution
                .calculate_helmholtz_energy_density(&density.mapv(Dual64::from), contributions)?,
        );

        for (c, wd) in functional_contributions.iter().zip(weighted_densities) {
            let nwd = wd.shape()[0];
            let ngrid = wd.len() / nwd;
            helmholtz_energy_density.push(
                c.calculate_helmholtz_energy_density(
                    temperature_dual,
                    wd.into_shape((nwd, ngrid)).unwrap().view(),
                )?
                .into_shape(density.raw_dim().remove_axis(Axis(0)))
                .unwrap(),
            );
        }
        Ok(helmholtz_energy_density
            .iter()
            .map(|v| v.mapv(|f| -(f * temperature_dual).eps[0]))
            .collect())
    }

    pub fn internal_energy_density<D>(
        &self,
        temperature: f64,
        density: &Array<f64, D::Larger>,
        external_potential: &Array<f64, D::Larger>,
        convolver: &Rc<dyn Convolver<Dual64, D>>,
        contributions: Contributions,
    ) -> EosResult<Array<f64, D>>
    where
        D: Dimension,
        D::Larger: Dimension<Smaller = D>,
    {
        let temperature_dual = Dual64::from(temperature).derive();
        let helmholtz_energy_density_dual = self.intrinsic_helmholtz_energy_density(
            temperature_dual,
            density,
            convolver,
            contributions,
        )?;
        let helmholtz_energy_density = helmholtz_energy_density_dual
            .mapv(|f| f.re - f.eps[0] * temperature)
            + (external_potential * density).sum_axis(Axis(0)) * temperature;
        Ok(helmholtz_energy_density)
    }

    #[allow(clippy::type_complexity)]
    pub fn functional_derivative<D>(
        &self,
        temperature: f64,
        density: &Array<f64, D::Larger>,
        convolver: &Rc<dyn Convolver<f64, D>>,
    ) -> EosResult<(Array<f64, D>, Array<f64, D::Larger>)>
    where
        D: Dimension,
        D::Larger: Dimension<Smaller = D>,
    {
        let weighted_densities = convolver.weighted_densities(density);
        let contributions = self.functional.contributions();
        let mut partial_derivatives = Vec::with_capacity(contributions.len());
        let mut helmholtz_energy_density = Array::zeros(density.raw_dim().remove_axis(Axis(0)));
        for (c, wd) in contributions.iter().zip(weighted_densities) {
            let nwd = wd.shape()[0];
            let ngrid = wd.len() / nwd;
            let mut phi = Array::zeros(density.raw_dim().remove_axis(Axis(0)));
            let mut pd = Array::zeros(wd.raw_dim());
            c.first_partial_derivatives(
                temperature,
                wd.into_shape((nwd, ngrid)).unwrap(),
                phi.view_mut().into_shape(ngrid).unwrap(),
                pd.view_mut().into_shape((nwd, ngrid)).unwrap(),
            )?;
            partial_derivatives.push(pd);
            helmholtz_energy_density += &phi;
        }
        Ok((
            helmholtz_energy_density,
            convolver.functional_derivative(&partial_derivatives),
        ))
    }

    // iSAFT correction to the functional derivative
    pub fn isaft_integrals<D>(
        &self,
        temperature: f64,
        functional_derivative: &Array<f64, D::Larger>,
        convolver: &Rc<dyn Convolver<f64, D>>,
    ) -> Array<f64, D::Larger>
    where
        D: Dimension,
        D::Larger: Dimension<Smaller = D>,
    {
        // calculate weight functions
        let bond_lengths = self.functional.bond_lengths(temperature).into_edge_type();
        let mut isaft_weight_functions = bond_lengths.map(
            |_, _| (),
            |_, &l| WeightFunction::new_scaled(arr1(&[l]), WeightFunctionShape::Delta),
        );
        for n in bond_lengths.node_indices() {
            for e in bond_lengths.edges(n) {
                isaft_weight_functions.add_edge(
                    e.target(),
                    e.source(),
                    WeightFunction::new_scaled(arr1(&[*e.weight()]), WeightFunctionShape::Delta),
                );
            }
        }

        let expdfdrho = functional_derivative.mapv(|x| (-x).exp());
        let mut i_graph: Graph<_, Option<Array<f64, D>>, Directed> =
            isaft_weight_functions.map(|_, _| (), |_, _| None);

        let bonds = i_graph.edge_count();
        let mut calc = 0;

        // go through the whole graph until every bond has been calculated
        while calc < bonds {
            let mut edge_id = None;
            let mut i1 = None;

            // find the first bond that can be calculated
            'nodes: for node in i_graph.node_indices() {
                for edge in i_graph.edges(node) {
                    // skip already calculated bonds
                    if edge.weight().is_some() {
                        continue;
                    }

                    // if all bonds from the neighboring segment are calculated calculate the bond
                    let edges = i_graph
                        .edges(edge.target())
                        .filter(|e| e.target().index() != node.index());
                    if edges.clone().all(|e| e.weight().is_some()) {
                        edge_id = Some(edge.id());
                        let i0 = edges.fold(
                            expdfdrho
                                .index_axis(Axis(0), edge.target().index())
                                .to_owned(),
                            |acc: Array<f64, D>, e| acc * e.weight().as_ref().unwrap(),
                        );
                        i1 = Some(
                            convolver.convolve(i0.clone(), &isaft_weight_functions[edge.id()]),
                        );
                        break 'nodes;
                    }
                }
            }
            if let Some(edge_id) = edge_id {
                i_graph[edge_id] = i1;
                calc += 1;
            } else {
                panic!("Cycle in molecular structure detected!")
            }
        }

        let mut i = Array::ones(functional_derivative.raw_dim());
        for node in i_graph.node_indices() {
            for edge in i_graph.edges(node) {
                i.index_axis_mut(Axis(0), node.index())
                    .mul_assign(edge.weight().as_ref().unwrap());
            }
        }

        i
    }
}
