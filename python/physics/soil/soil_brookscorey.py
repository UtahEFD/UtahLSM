from .soil import SOIL

class BROOKSCOREY(SOIL):

    # class initialization
    def __init__(self,inputLSM):
        
        print("[UtahLSM: SOIL] \t--- running with the Brooks-Corey scheme")
        # initialize parent class
        super().__init__(inputLSM)