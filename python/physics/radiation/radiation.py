class Radiation(object):
    
    def __init__(self,inputLSM):
    
        self.input      = inputLSM
        self.latitude   = inputLSM.latitude
        self.longitude  = inputLSM.longitude
        self.albedo     = inputLSM.albedo
        self.emissivity = inputLSM.emissivity
        
    @staticmethod
    def get_model(key,inputLSM):
        if key == 1:
            from .rad_basic import RadBasic
            return RadBasic(inputLSM)