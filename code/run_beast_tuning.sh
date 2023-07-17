for steps in [4,8,16,32,64]; do
	for rates in [0.2,0.3,0.4,0.5,0.6,0.7,0.8]; do
		java -jar -Djava.library.path=/home/andrew/lib /home/andrew/beast-mcmc/build/dist/beast.jar -seed 666 -overwrite /home/andrew/mcmc_handbook/code/tuning$steps$rates.xml &
	done
done