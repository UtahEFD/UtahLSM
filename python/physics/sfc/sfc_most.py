from .sfc import SFC

class MOST(SFC):

    # class initialization
    def __init__(self,inputLSM):
        
        print("[UtahLSM: SL] \t\t --- running with the MOST scheme")
        # initialize parent class
        super().__init__(inputLSM)