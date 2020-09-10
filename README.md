# CTS-Ext
# A full-scale agent-based model of Lombardy COVID-19 outbreak to explore social networks connectivity

***Author***: **Giuseppe Giacopelli** [ORCID: 0000-0001-6801-8811]

***Abstarct***: COVID-19 outbreak is an awful event. However it gives to the scientists the possibility to test theories about epidemic. The aim of this contribution is to propose a individual-based model of Lombardy COVID-19 outbreak at full-scale, where full-scale means that will be simulated all the 10 millions inhabitant population of Lombardy person by person, in a commercial computer. All this to test the impact of our daily actions in epidemic and in the end have an insight on social networks connectivity.

***Usage***: This project is MATLAB based and has three main scripts: The first one (_lombardia.m_) runs the epidemic simulation, the scond one (_deg.m_) computes the degree distribution of a small group of people, the third one (_test_acc.m_) estimates the prcision of collision detection algorithm. You have just to launch the script and wait for the results. The time scale for the scripts is about 30 minutes for a machine with 64GB of RAM and an AMD 3900X 12-cores processor. Inside the the script the variable _mode_ select the mode of the simulation. Main values for mode are: _st_ is the standard mode, _ld_ the lock down mode and _2m_ the social distancing mode. Uncomment the value that you prefer.
