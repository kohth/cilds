function numWorkers = StartParPool(numParRuns)
distcomp.feature( 'LocalUseMpiexec', false );
currParPool = gcp('NoCreate');

numWorkers = min(feature('NumCores'), numParRuns);

if isempty(currParPool)
	parpool(numWorkers);
elseif currParPool.NumWorkers < numWorkers
	delete(currParPool);
	parpool(numWorkers);
end

end
