from .radiation import RAD

class BASIC(RAD):

    # class initialization
    def __init__(self,inputLSM):
        
        print("[UtahLSM: RAD] \t\t --- running with the basic model")
        # initialize parent class
        super().__init__(inputLSM)