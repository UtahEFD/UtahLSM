class Surface(object):
    
    def __init__(self):pass
        
    @staticmethod
    def get_model(key):
        if key == 1:
            from .sfc_most import SurfaceMOST
            return SurfaceMOST()