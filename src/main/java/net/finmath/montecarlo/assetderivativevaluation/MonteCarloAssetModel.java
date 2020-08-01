/*
 * (c) Copyright Christian P. Fries, Germany. Contact: email@christian-fries.de.
 *
 * Created on 20.01.2004
 */
package net.finmath.montecarlo.assetderivativevaluation;

import java.time.LocalDateTime;
import java.util.Map;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.model.AbstractProcessModel;
import net.finmath.montecarlo.model.ProcessModel;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.stochastic.RandomVariable;
import net.finmath.time.TimeDiscretization;

/**
 * This class glues together an <code>AbstractProcessModel</code> and a Monte-Carlo implementation of a <code>MonteCarloProcessFromProcessModel</code>
 * and implements <code>AssetModelMonteCarloSimulationModel</code>.
 *
 * The model is specified via the object implementing <code>ProcessModel</code>.
 *
 * @author Christian Fries
 * @see net.finmath.montecarlo.process.MonteCarloProcess The interface for numerical schemes.
 * @see net.finmath.montecarlo.model.ProcessModel The interface for models provinding parameters to numerical schemes.
 * @version 1.0
 */
public class MonteCarloAssetModel implements AssetModelMonteCarloSimulationModel {

	private final ProcessModel model;

	/**
	 * Create a Monte-Carlo simulation using given process discretization scheme.
	 *
	 * @param model The model to be used.
	 * @param process The numerical scheme to be used.
	 */
	public MonteCarloAssetModel(
			final ProcessModel model,
			final MonteCarloProcess process) {
		super();

		this.model = model;

		// Link model and process for delegation
		process.setModel(model);
		model.setProcess(process);
	}

	@Override
	public RandomVariable getAssetValue(final double time, final int assetIndex) throws CalculationException {
		return getAssetValue(getTimeIndex(time), assetIndex);
	}

	@Override
	public RandomVariable getAssetValue(final int timeIndex, final int assetIndex) throws CalculationException {
		return model.getProcess().getProcessValue(timeIndex, assetIndex);
	}

	@Override
	public RandomVariable getNumeraire(final int timeIndex) throws CalculationException {
		final double time = getTime(timeIndex);

		return model.getNumeraire(time);
	}

	@Override
	public RandomVariable getNumeraire(final double time) throws CalculationException {
		return model.getNumeraire(time);
	}

	@Override
	public RandomVariable getMonteCarloWeights(final double time) throws CalculationException {
		return getMonteCarloWeights(getTimeIndex(time));
	}

	@Override
	public int getNumberOfAssets() {
		return 1;
	}

	@Override
	public AssetModelMonteCarloSimulationModel getCloneWithModifiedData(final Map<String, Object> dataModified) throws CalculationException {
		final MonteCarloProcess process = model.getProcess();

		final ProcessModel		newModel	= model.getCloneWithModifiedData(dataModified);

		MonteCarloProcess	newProcess;
		try {
			newProcess = process.getCloneWithModifiedData(dataModified);
		}
		catch(final UnsupportedOperationException e) {
			newProcess = process;
		}

		// In the case where the model has changed we need a new process anyway
		if(newModel != model && newProcess == process) {
			newProcess = process.clone();
		}

		return new MonteCarloAssetModel(newModel, newProcess);
	}

	@Override
	public AssetModelMonteCarloSimulationModel getCloneWithModifiedSeed(final int seed) {
		throw new UnsupportedOperationException("Method not implemented");
	}

	@Override
	public int getNumberOfPaths() {
		return model.getProcess().getNumberOfPaths();
	}

	@Override
	public LocalDateTime getReferenceDate() {
		return model.getReferenceDate();
	}

	@Override
	public TimeDiscretization getTimeDiscretization() {
		return model.getProcess().getTimeDiscretization();
	}

	/* (non-Javadoc)
	 * @see net.finmath.montecarlo.MonteCarloSimulationModel#getTime(int)
	 */
	@Override
	public double getTime(final int timeIndex) {
		return model.getProcess().getTime(timeIndex);
	}

	/* (non-Javadoc)
	 * @see net.finmath.montecarlo.MonteCarloSimulationModel#getTimeIndex(double)
	 */
	@Override
	public int getTimeIndex(final double time) {
		return model.getProcess().getTimeIndex(time);
	}

	@Override
	public RandomVariable getRandomVariableForConstant(final double value) {
		return model.getRandomVariableForConstant(value);
	}

	@Override
	public RandomVariable getMonteCarloWeights(final int timeIndex) throws CalculationException {
		return model.getProcess().getMonteCarloWeights(timeIndex);
	}

	/**
	 * Returns the {@link AbstractProcessModel} used for this Monte-Carlo simulation.
	 *
	 * @return the model
	 */
	public ProcessModel getModel() {
		return model;
	}

	@Override
	public String toString() {
		return this.getClass().getSimpleName() + " [model=" + model + "]";
	}
}
