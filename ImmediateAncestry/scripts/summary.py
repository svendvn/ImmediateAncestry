class Summary(object):
       
    def __init__(self, name, output='double'):
        self.name=name
        self.output=output
    
    def __call__(self, **kwargs):
        pass
    
    def pretty_print(self, output):
        return self.name+'= '+str(output)
    

class basic_config_summary(object):
    
    def __init__(self, function, name, output='double'):
        super(basic_config_summary, self).__init__(name, output=output)
        self.function=function
        
    def __call__(self, **kwargs):
        tree=kwargs['tree']
        return self.function(tree)

    def summary_of_phylogeny(self, tree):
        return self.function(tree)