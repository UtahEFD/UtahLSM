class RAD(object):
    
    def __init__(self,inputLSM):
    
        self.input = inputLSM
        
    @staticmethod
    def get_model(key,inputLSM):
        if key == 1:
            from .rad_basic import BASIC
            return BASIC(inputLSM)