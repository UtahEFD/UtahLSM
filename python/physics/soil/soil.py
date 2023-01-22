class SOIL(object):
    
    def __init__(self,inputLSM):
    
        self.input = inputLSM
        
    @staticmethod
    def get_model(key,inputLSM):
        if key == 1:
            from .soil_brookscorey import BROOKSCOREY
            return BROOKSCOREY(inputLSM)
        if key == 2:
            from .soil_campbell import CAMPBELL
            return CAMPBELL(inputLSM)
        if key == 3:
            from .soil_vangenuchten import VANGENUCHTEN
            return VANGENUCHTEN(inputLSM)