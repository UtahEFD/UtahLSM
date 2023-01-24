class Radiation(object):
    
    def __init__(self, input):
    
        self.latitude   = input.latitude
        self.longitude  = input.longitude
        self.albedo     = input.albedo
        self.emissivity = input.emissivity
        
    @staticmethod
    def get_model(key, input):
        if key == 1:
            from .rad_basic import RadBasic
            return RadBasic(input)