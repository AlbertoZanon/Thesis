/*
 * (c) Copyright Christian P. Fries, Germany. Contact: email@christian-fries.de.
 *
 * Created on 22.12.2016
 */

package net.finmath.montecarlo.interestrate.models.covariance;

import net.finmath.montecarlo.interestrate.TermStructureModel;
import net.finmath.montecarlo.interestrate.models.LIBORMarketModelWithTenorRefinement;
import net.finmath.stochastic.RandomVariable;
import net.finmath.time.TimeDiscretization;

/**
 * @author Christian Fries
 *
 * @version 1.0
 */
public class TermStructCovarianceModelFromLIBORCovarianceModel implements TermStructureFactorLoadingsModelInterface {

	private final AbstractLIBORCovarianceModelParametric covarianceModel;

	/**
	 * Create a term structure covariance model model implementing TermStructureCovarianceModelInterface
	 * using a given model implementing AbstractLIBORCovarianceModelParametric.
	 *
	 * @param covarianceModel The model implementing AbstractLIBORCovarianceModelParametric.
	 */
	public TermStructCovarianceModelFromLIBORCovarianceModel(final AbstractLIBORCovarianceModelParametric covarianceModel) {
		this.covarianceModel = covarianceModel;
	}

	@Override
	public RandomVariable[] getFactorLoading(final double time, final double periodStart, final double periodEnd, final TimeDiscretization periodDiscretization, final RandomVariable[] realizationAtTimeIndex, final TermStructureModel model) {
		final TimeDiscretization liborPeriodDiscretization = covarianceModel.getLiborPeriodDiscretization();

		// Cache is really needed.
		final RandomVariable[] liborAtTimeIndex = new RandomVariable[liborPeriodDiscretization.getNumberOfTimeSteps()];
		for(int componentIndex=0; componentIndex<liborPeriodDiscretization.getNumberOfTimeSteps(); componentIndex++) {
			if(liborPeriodDiscretization.getTime(componentIndex) < time) {
				liborAtTimeIndex[componentIndex] = null;
			}
			else {
				// in realtà anche questo metodo ritorna NULL. Da chiedere
				liborAtTimeIndex[componentIndex] = ((LIBORMarketModelWithTenorRefinement)model).getLIBORForStateVariable(periodDiscretization, realizationAtTimeIndex, liborPeriodDiscretization.getTime(componentIndex), liborPeriodDiscretization.getTime(componentIndex+1));
			}
		}

		final int periodStartIndex = liborPeriodDiscretization.getTimeIndex(periodStart);
		final int periodEndIndex = liborPeriodDiscretization.getTimeIndex(periodEnd);
//Quindi anche qua gli passa sempre liborAtTimeIndex=null --> se deve fare blended covariance?		
		final RandomVariable[] factorLoadings = covarianceModel.getFactorLoading(time, periodStartIndex, liborAtTimeIndex);
		if(periodEndIndex > periodStartIndex+1) {
			// Need to sum factor loadings
			for(int factorIndex = 0; factorIndex<factorLoadings.length; factorIndex++) {
				//lambda_i=lambda_i*deltaT
				factorLoadings[factorIndex] = factorLoadings[factorIndex].mult(liborPeriodDiscretization.getTimeStep(periodStartIndex));
			}
			
			for(int periodIndex = periodStartIndex+1; periodIndex<periodEndIndex; periodIndex++) {
				final RandomVariable[] factorLoadingsForPeriod = covarianceModel.getFactorLoading(time, periodStartIndex, liborAtTimeIndex);
				final double periodLength = liborPeriodDiscretization.getTimeStep(periodIndex);
				for(int factorIndex = 0; factorIndex<factorLoadings.length; factorIndex++) {
					// lambda_i = lambda_i(periodIndex) + lambda_i(periodIndex)*deltaT + . . . (next loop) . . . + lambda_i(periodIndex + 1)*deltaT
					factorLoadings[factorIndex] = factorLoadings[factorIndex].addProduct(factorLoadingsForPeriod[factorIndex], periodLength);
				}
			}

			for(int factorIndex = 0; factorIndex<factorLoadings.length; factorIndex++) {
				//lambda_i = lambda_i / (periodEnd-periodStart)
				factorLoadings[factorIndex] = factorLoadings[factorIndex].div(periodEnd-periodStart);
			}
		}
		return factorLoadings;
	}

	@Override
	public int getNumberOfFactors() {
		return covarianceModel.getNumberOfFactors();
	}
}
