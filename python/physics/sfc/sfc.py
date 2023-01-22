class SFC(object):
    
    def __init__(self,inputLSM):
    
        self.input = inputLSM
        
    @staticmethod
    def get_model(key,inputLSM):
        if key == 1:
            from .sfc_most import MOST
            return MOST(inputLSM)