from .soil import SOIL

class CAMPBELL(SOIL):

    # class initialization
    def __init__(self,inputLSM):
        
        print("[UtahLSM: SOIL] \t--- running with the Campbell scheme")
        # initialize parent class
        super().__init__(inputLSM)