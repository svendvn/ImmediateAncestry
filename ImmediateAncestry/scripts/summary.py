class Summary(object):
       
    def __init__(self, name, output='double'):
        self.name=name
        self.output=output
    
    def __call__(self, **kwargs):
        pass
    
    def pretty_print(self, output):
        return self.name+'= '+str(output)
    

class basic_config_summary(Summary):
    
    def __init__(self, function, name, output='double'):
        super(basic_config_summary, self).__init__(name, output=output)
        self.function=function
        
    def __call__(self, **kwargs):
        x=kwargs['x']
        return self.function(x)

    def summary_of_x(self, z):
        return self.function(z)
    
class basic_summary(Summary):
    
    def __init__(self, name, output='double'):
        super(basic_summary, self).__init__(name, output=output)
        
    def __call__(self, **kwargs):
        x=kwargs[self.name]
        return x

    