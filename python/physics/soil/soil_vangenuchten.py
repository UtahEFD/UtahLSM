from .soil import SOIL

class VANGENUCHTEN(SOIL):

    # class initialization
    def __init__(self,inputLSM):
        
        print("[UtahLSM: SOIL] \t--- running with the Van Genuchten scheme")
        # initialize parent class
        super().__init__(inputLSM)