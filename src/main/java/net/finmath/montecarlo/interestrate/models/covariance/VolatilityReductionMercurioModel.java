/*
 * (c) Copyright Christian P. Fries, Germany. Contact: email@christian-fries.de.
 *
 * Created on 26.05.2013
 */
package net.finmath.montecarlo.interestrate.models.covariance;

import java.util.Map;

import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.montecarlo.RandomVariableFactory;
import net.finmath.stochastic.RandomVariable;
import net.finmath.stochastic.Scalar;

/**
 * Volatility reduction model build on top of a standard covariance model.
 *
 * @param covarianceModel The given covariance model specifying the factor loadings <i>F</i>.
 * @param isCalibrateable Note that the covariance model may have its own parameter calibration settings.
 */
public class VolatilityReductionMercurioModel extends AbstractLIBORCovarianceModelParametric {

	private static final long serialVersionUID = 4522227972747028512L;
	private final AbstractLIBORCovarianceModelParametric covarianceModel;

	/**
	 * Volatility reduction model build on top of a standard covariance model.
	 *
	 * @param covarianceModel The given covariance model specifying the factor loadings <i>F</i>.
	 */
	public VolatilityReductionMercurioModel(final AbstractLIBORCovarianceModelParametric covarianceModel) {
		super(covarianceModel.getTimeDiscretization(), covarianceModel.getLiborPeriodDiscretization(), covarianceModel.getNumberOfFactors());
		this.covarianceModel	= covarianceModel;
	}


	@Override
	public Object clone() {
		return new VolatilityReductionMercurioModel((AbstractLIBORCovarianceModelParametric) covarianceModel.clone());
	}

	/**
	 * Returns the base covariance model, i.e., the model providing the factor loading <i>F</i>
	 * such that this model's <i>i</i>-th factor loading is
	 * <i>(a L<sub>i,0</sub> + (1-a)L<sub>i</sub>(t)) F<sub>i</sub>(t)</i>
	 * where <i>a</i> is the displacement and <i>L<sub>i</sub></i> is
	 * the realization of the <i>i</i>-th component of the stochastic process and
	 * <i>F<sub>i</sub></i> is the factor loading loading from the given covariance model.
	 *
	 * @return The base covariance model.
	 */
	public AbstractLIBORCovarianceModelParametric getBaseCovarianceModel() {
		return covarianceModel;
	}

	@Override
	public RandomVariable[] getParameter() {
		return covarianceModel.getParameter();
	}

	@Override
	public double[] getParameterAsDouble() {
		final RandomVariable[] parameters = getParameter();
		final double[] parametersAsDouble = new double[parameters.length];
		for(int i=0; i<parameters.length; i++) {
			parametersAsDouble[i] = parameters[i].doubleValue();
		}
		return parametersAsDouble;
	}

	@Override
	public AbstractLIBORCovarianceModelParametric getCloneWithModifiedParameters(final RandomVariable[] parameters) {
		if(parameters == null || parameters.length == 0) {
			return this;
		}
			return new VolatilityReductionMercurioModel(covarianceModel.getCloneWithModifiedParameters(parameters));
	}

	@Override
	public AbstractLIBORCovarianceModelParametric getCloneWithModifiedParameters(final double[] parameters) {
		return getCloneWithModifiedParameters(Scalar.arrayOf(parameters));
	}
	

	@Override
	public RandomVariable[] getFactorLoading(final int timeIndex, final int component, final RandomVariable[] realizationAtTimeIndex) {
		final RandomVariable[] factorLoading = covarianceModel.getFactorLoading(timeIndex, component, realizationAtTimeIndex);

		if(realizationAtTimeIndex != null && realizationAtTimeIndex[component] != null) {
			final double simulationTime =  getTimeDiscretization().getTime(timeIndex);
			final double liborEndTime = getLiborPeriodDiscretization().getTime(component+1);
			final double liborStartTime = getLiborPeriodDiscretization().getTime(component);
			
/*
 *  		apply volatility reduction function on the factorLoading of covarianceModel.
 *   		volatility reduction function g(t)=min(\frac{(T_{i+1}-t)^+}{T_{j+1}-T_j},1) 
 *   		the reduction thus apply only on the accrual period [T_i,T_{i+1}], otherwise = 1.
 */			
			final double Volatilityreduction = Math.min(Math.max(liborEndTime-simulationTime, 0)/(liborEndTime-liborStartTime), 1);		
			for (int factorIndex = 0; factorIndex < factorLoading.length; factorIndex++) {
				factorLoading[factorIndex] = factorLoading[factorIndex].mult(Volatilityreduction);
			}
		}

		return factorLoading;
	}

	@Override
	public RandomVariable getFactorLoadingPseudoInverse(final int timeIndex, final int component, final int factor, final RandomVariable[] realizationAtTimeIndex) {
		throw new UnsupportedOperationException();
	}

	@Override
	public AbstractLIBORCovarianceModelParametric getCloneWithModifiedData(final Map<String, Object> dataModified)
			throws CalculationException {

		AbstractLIBORCovarianceModelParametric covarianceModel = this.covarianceModel;

		if(dataModified != null) {
			if (!dataModified.containsKey("covarianceModel")) {
				covarianceModel = covarianceModel.getCloneWithModifiedData(dataModified);
			}
			// Explicitly passed covarianceModel has priority
			covarianceModel = (AbstractLIBORCovarianceModelParametric)dataModified.getOrDefault("covarianceModel", covarianceModel);

		}
		final AbstractLIBORCovarianceModelParametric newModel = new VolatilityReductionMercurioModel(covarianceModel);
		return newModel;
	}
}
