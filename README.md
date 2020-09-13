# CTS-Ext
# A full-scale agent-based model of Lombardy COVID-19 outbreak to explore social networks connectivity

***Author***: **Giuseppe Giacopelli** [ORCID: 0000-0001-6801-8811]

***Abstarct***: COVID-19 outbreak is an awful event. However it gives to the scientists the possibility to test theories about epidemic. The aim of this contribution is to propose a individual-based model of Lombardy COVID-19 outbreak at full-scale, where full-scale means that will be simulated all the 10 millions inhabitant population of Lombardy person by person, in a commercial computer. All this to test the impact of our daily actions in epidemic and in the end have an insight on social networks connectivity.

***Usage***: This project is MATLAB based and has three main scripts: The first one (_lombardia.m_) runs the epidemic simulation, the scond one (_deg.m_) computes the degree distribution of a small group of people, the third one (_test_acc.m_) estimates the prcision of collision detection algorithm. You have just to launch the script and wait for the results. The time scale for the scripts is about 30 minutes for a machine with 64GB of RAM and an AMD 3900X 12-cores processor. Inside each script the variable _mode_ selects the mode of simulation. Main values for _mode_ are: _st_ for standard mode, _ld_ for lock down mode and _2m_ for social distancing mode. Uncomment the value that you prefer.

***Acknowledgement***: The author received no financial support for the research of this article and all the simulations have been performed on the personal computer of the author. For these reasons doesn't exist any conflict of interest in publishing the source code of the project here. However author thanks the department of Mathematics and Informatics of the University of Palermo and the Institute of Biophysics of National Research Council (CNR) for funding author's PhD and for being a constant guide for his growth.
