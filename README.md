# RNaiveBayesianModeler

R source file to import raw data files, and generate a Bayesian model from the data.

FUNCTIONS:

	binCreator(data.frame df, array subjMat)
		Bin creator accepts a data.frame of raw data.  It uses an array of
		subject IDs (subjMat) to create bins for each electrode for each subject
		in subjMat.

		Returns data.frame of bins of electrode potentials.

	binFrequency
