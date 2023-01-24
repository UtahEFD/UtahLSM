class Radiation(object):
    
    def __init__(self, latitude, longitude, albedo, emissivity):
    
        self.latitude   = latitude
        self.longitude  = longitude
        self.albedo     = albedo
        self.emissivity = emissivity
        
    @staticmethod
    def get_model(key,latitude, longitude, albedo, emissivity):
        if key == 1:
            from .rad_basic import RadBasic
            return RadBasic(latitude, longitude, albedo, emissivity)