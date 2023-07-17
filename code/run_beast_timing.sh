for steps in {4,8,16,32,64}; do
	time java -jar -Djava.library.path=/home/andrew/lib /home/andrew/beast-mcmc/build/dist/beast.jar -seed 666 -overwrite /home/andrew/mcmc_handbook/code/timing$steps.xml 
done