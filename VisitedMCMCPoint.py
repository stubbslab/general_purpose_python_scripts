#Defines the VisitedDsphMCMCPoint class.
#This class is used to keep track of each point visited in the parameter space during each MCMC chain.
#It is composed of a few pieces :
#  the parameter storer, which holds all of the parameters
#  the log likelihood of the visited point (based on the normalized surface brightness profile, given the parameters used)
#  the number of times the algorithm stayed at the point (starting at 1, and increasing by 1 every time the algorithm doesn't move from a point)
#Then when, the MCMC series is run, the output consists of an array of instances of this class.


class VisitedMCMCPoint:

    def __init__(self, params, likelihood, n_visits):
        self.parameters= params
        self.likelihood = likelihood
        self.n_visits = n_visits
